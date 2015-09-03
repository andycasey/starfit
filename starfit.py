import numpy as np
import cPickle as pickle
import os
from collections import Counter,OrderedDict

import requests
import xlrd
from requests_toolbelt import MultipartEncoder


OVERWRITE = False
DEFAULT_UNCERTAINTY = 0.2
FILENAME = "data/literature.xlsx"
RESULTS_FOLDER = "results-Zmax_30"


STARFIT_URL = "http://starfit.org/py/starfitonline.py"
STARFIT_DEFAULT_PAYLOAD = OrderedDict([
    # Default database: znuc2012.S4.star.el.y.stardb.gz
    ("algorithm", "single"),
    ("time_limit", "5"), 
    ("pop_size", "200"),
    ("gene_size", "2"),     
    ("fixed", "0"),
    ("stardata", None),
    ("database", "znuc2012.S4.star.el.y.stardb.gz"),
    ("z_lolim", "Sc, Cu"),
    ("combine", "0"),         
    ("z_max", "30"),
    ("z_exclude", "Li"), #Online defaults: Li, Cr, Zn
    ("email", ""),
    ("plotformat", "svg"),
])



def parse_starfit_result(content):
    
    def get_best_model(string):
        _ = string.split(": mass = ")
        model_id = int(_[0].split("\n")[-1].strip())

        _ = _[1].split("\n")[0]
        params = _.split(", ")
        params[1:] = [each.split(" = ")[1] for each in params[1:]]
        mass, energy, log_mixing, remnant = map(float, params)

        print("Returned model id {0} with mass {1:.1f}, energy {2:.2f} and remnant mass {3:.2f}".format(
            model_id, mass, energy, remnant))
        
        return {
            "model_id": model_id,
            "mass": mass,
            "energy": energy,
            "remnant": remnant,
        }

    def get_mean_squared_residual(string):
        return {
            "mean_squared_residual": float(string.split("Mean squared residual: </b>")[1].split("\n")[0].strip())
        }


    result = {}
    result.update(get_best_model(content))
    result.update(get_mean_squared_residual(content))
    return result


def atomic_number(element):
    """
    Return the atomic number of a given element.

    :param element:
        The short-hand notation for the element (e.g., Fe).

    :type element:
        str

    :returns:
        The atomic number for a given element.

    :rtype:
        int
    """
    
    if not isinstance(element, (unicode, str)):
        raise TypeError("element must be represented by a string-type")

    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    
    if element not in periodic_table:
        return ValueError("element '{0}' is not known".format(element))

    return periodic_table.index(element) + 1



def starfit_format(star_name, perturb=False, adjust_C=False):

    if adjust_C:
        raise NotImplementedError


    hdr_row_index = 18

    sheet = data.sheet_by_name(star_name)

    literature_source = sheet.row(3)[1].value

    ra = sheet.row(4)[1].value
    dec = sheet.row(5)[1].value

    # Get teff and logg
    teff = sheet.row(7)[1].value
    logg = sheet.row(8)[1].value

    use_abundances = sheet.row(12)[1].value
    solar_abundances = sheet.row(14)[1].value
    comment = sheet.row(17)[1].value

    data_format = {
    "log_eps": 1,
    "[X/H": 6,
    "[X/H]": 6
    }.get(use_abundances, None)
    assert data_format is not None

    # Now extract all abundances.
    uncertainty_column_index = 4
    _ = [_.value for _ in sheet.row(hdr_row_index)]
    assert use_abundances in _
    abundance_column_index = _.index(use_abundances)

    C_measurement_exists = False
    N_measurement_exists = False
    O_measurement_exists = False
    abundance_info = {}
    for row_index in range(hdr_row_index + 1, sheet.nrows):

        row = sheet.row(row_index)

        element = row[0].value
        abundance = str(row[abundance_column_index].value)

        # Check if it's an upper limit, and return the value.
        is_upper_limit = "<" in abundance
        abundance = abundance.replace("<", "").strip()

        # If no information is given, just continue.
        if len(abundance) == 0:
            continue

        if element == "C" and not is_upper_limit:
            C_measurement_exists = True

        if element == "N" and not is_upper_limit:
            N_measurement_exists = True

        if element == "O" and not is_upper_limit:
            O_measurement_exists = True

        uncertainty = str(row[uncertainty_column_index].value)

        if is_upper_limit:
            uncertainty = "-0.20" # As per the SM0313-6708 example.

        if len(uncertainty.strip()) == 0:
            print("Assuming uncertainty of {0:.1f} for {1} abundance".format(
                DEFAULT_UNCERTAINTY, element))
            uncertainty = "+0.20"

        if perturb:
            if not is_upper_limit:
                abundance, uncertainty = map(float, (abundance, uncertainty))
                new_abundance = np.random.normal(abundance, uncertainty)
                print("Perturbing measured {0} abundance from {1:.2f} +/- {2:.2f} to be {3:.2f}".format(element, abundance, uncertainty, new_abundance))
                abundance = new_abundance

            else:
                # [X/H] drawn from uniform (-6, X)
                # Need abundances to translate, since we have a limit in log(eps)
                # and we want a value in [X/H] = log_eps - log_eps_solar
                # log_eps = [X/H] + log_eps_solar
                __solar = {
                    "Li": 1.05,
                    "Be": 1.38,
                    "B": 2.70,
                    "C": 8.43,
                    "N": 7.83,
                    "O": 8.69,
                    "F": 4.56,
                    "Ne": 7.93,
                    "Na": 6.24,
                    "Mg": 7.60,
                    "Al": 6.45,
                    "Si": 7.51,
                    "P": 5.41,
                    "S": 7.12,
                    "Cl": 5.50,
                    "Ar": 6.40,
                    "K": 5.03,
                    "Ca": 6.34,
                    "Sc": 3.15,
                    "Ti": 4.95,
                    "V": 3.93,
                    "Cr": 5.64,
                    "Mn": 5.43,
                    "Fe": 7.50,
                    "Co": 4.99,
                    "Ni": 6.22,
                    "Cu": 4.19,
                    "Zn": 4.56,
                    "Ga": 3.04,
                    "Ge": 3.65,
                    "Kr": 3.25,
                    "Rb": 2.52,
                    "Sr": 2.87,
                    "Y": 2.21,
                    "Zr": 2.58,
                    "Nb": 1.46,
                    "Mo": 1.88,
                    "Ru": 1.75,
                    "Rh": 0.91,
                    "Pd": 1.57,
                    "Ag": 0.94,
                    "In": 0.80,
                    "Sn": 2.04,
                    "Ba": 2.18,
                    "La": 1.10,
                    "Ce": 1.58,
                    "Pr": 0.72,
                    "Nd": 1.42,
                    "Sm": 0.96,
                    "Eu": 0.52,
                    "Gd": 1.07,
                    "Tb": 0.30,
                    "Dy": 1.10,
                    "Ho": 0.48,
                    "Er": 0.92,
                    "Tm": 0.10,
                    "Yb": 0.84,
                    "Lu": 0.10,
                    "Hf": 0.85,
                    "W": 0.85,
                    "Os": 1.40,
                    "Ir": 1.38,
                    "Au": 0.92,
                    "Tl": 0.90,
                    "Pb": 1.75,
                    "Th": 0.02
                }.get(element.split(" ")[0])
                x_h_abundance = float(abundance) - __solar
                new_abundance = np.random.uniform(-6, x_h_abundance)
                abundance = new_abundance + __solar
                
                uncertainty = uncertainty.replace("-", "+")

                print("Perturbing upper limit of {0} abundance to be drawn from [{0}/H] ~ Uniform(-6, {1:.2f}) to be {2:.2f} (log eps = {3:.2f})".format(element,
                    x_h_abundance, new_abundance, abundance))


        abundance_info[element] = (abundance, uncertainty)

    
    combine = 0
    if C_measurement_exists and N_measurement_exists:
        combine = 1
        if O_measurement_exists:
            combine = 2

    # Adjust C
    #if adjust_C:


    # Clean up duplicate elements.
    elements = [_.split(" ")[0] for _ in abundance_info.keys()]
    duplicate_elements = [k for k, v in Counter(elements).items() if v > 1]

    for element in duplicate_elements:
        del abundance_info["{} 1".format(element)]
        print("Taking {0} 2 as abundance over {0} 1".format(element))

    # Clean up the names
    z_max = 0
    abundances = {}
    upper_limits = []
    for k, v in abundance_info.items():
        k = k.split(" ")[0]
        if k.strip() == "": continue

        abundances[k] = v
        z_max = max([z_max, atomic_number(k)])
        if 0 > float(v[1]):
            upper_limits.append(k)

    assert isinstance(z_max, int)

    # Format the abundances.
    abundances_string = ""
    for k, v in abundances.items():
        abundances_string += "{0:5s} {1:6.2f} {2:6.2f}\n".format(k, *map(float, v))

    string = """
{version}
{star_name}
{literature_source}
{comment}
{data_format}
{abundance_offset}
{num_abundances}
{all_abundances}
{solar_abundance_source}

""".format(version=100002,
        star_name=star_name.replace("_", "-"),
        literature_source=literature_source,
        comment=comment,
        data_format=data_format, # assumed log_eps
        abundance_offset=0,
        num_abundances=len(abundances),
        all_abundances=abundances_string.rstrip(),
        solar_abundance_source="As09" # DEFAULT.
        )[1:]

    return (string, abundances.keys(), upper_limits, z_max, combine)



def starfit(contents, z_lolim=None, z_exclude=None, z_max=None, combine_CNO=0,
    **kwargs):
    """
    contents should be a filename or string.
    """

    if os.path.exists(contents):
        filename = os.path.basename(contents)
        stardata = open(filename, "rb")

    else:
        filename = "dat.dat"
        stardata = contents

    payload = STARFIT_DEFAULT_PAYLOAD.copy()
    payload.update({
        "z_lolim": payload["z_lolim"] if z_lolim is None else (", ".join(z_lolim) if isinstance(z_lolim, (list, tuple)) else z_lolim),
        "z_exclude": payload["z_exclude"] if z_exclude is None else (", ".join(z_exclude) if isinstance(z_exclude, (list, tuple)) else z_exclude),
        "z_max": str(z_max) if z_max is not None else payload["z_max"],
        "combine": str(combine_CNO),
        "stardata": (filename, stardata, "application/octet-stream")
    })

    data = MultipartEncoder(payload,
        boundary="----WebKitFormBoundarybmHMTFnJke5z2d0N")

    r = requests.post(STARFIT_URL, headers={
            "Host": "starfit.org",
            "Connection": "keep-alive",
            "Cache-Control": "max-age=0",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Origin": "http://starfit.org",
            "Upgrade-Insecure-Requests": 1,
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/44.0.2403.130 Safari/537.36",
            "DNT": 1,
            "Referer": "http://starfit.org/", 
            "Accept-Encoding": "gzip, deflate",
            "Accept-Language": "en-US,en;q=0.8",
            "Content-Type": data.content_type

        }, data=data.to_string().replace("\r", ""))


    result = parse_starfit_result(r.content)
    pickleable_payload = payload.copy()
    pickleable_payload["stardata"] = filename

    return (result, pickleable_payload, r.content)




if __name__ == "__main__":

    N_PERTURBATIONS = 30
    OVERWRITE = False

    data = xlrd.open_workbook(FILENAME)

    stars = data.sheet_names()[2:]

    def clean_star_name(name):
        return name.replace("$", "").replace("~", "").strip()

    """
    placco_stars = [
        ("CD-38~245", 21.5, 0.3),
        ("SDSS~J2209-0028", 27.0, 0.3),
        ("HE~2139-5432", 28.0, 0.6),
        ("G77-61", 27.0, 0.3),
        ("CS~30336-049", 21.5, 0.3),
        ("HE~0057-5959", 27.0, 0.3),
        ("HE~2239-5019", 15.0, 10.0),
        ("HE~1310-0536", 10.9, 0.3),
        ("SDSS~J1204$+$1201", 10.6, 0.9),
        ("HE~0233-0343", 11.9, 0.3),
        ("SDSS~J1742$+$2531", 21.5, 0.3),
        ("SDSS~J1029$+$1729", 10.6, 0.9),
        ("SDSS~J1313$+$0019", 27.0, 0.3),
        ("SDSS~J1035$+$0641", 23.0, 0.6),
        ("HE~0107-5240", 20.5, 0.6),
        ("HE~1327-2326_SUBG", 21.5, 0.3),
        ("SMSS~J0313-6708", 41.0, 1.2),
    ]
    """

    for star in stars:
        print("Star {}".format(star))

        prefix = clean_star_name(star)
        if os.path.exists("{0}/{1}.pkl".format(RESULTS_FOLDER, prefix)) and not OVERWRITE:
            print("Skipping {0} because {1}/{2}.pkl exists..".format(
                star, RESULTS_FOLDER, prefix))
            continue

        try:
            input_contents, elements, upper_limits, z_max, combine_CNO = starfit_format(star)

        except AssertionError:
            print("Skipping {0} because log_eps not specified for abundances".format(star))
            continue

        except:
            raise

        result, payload, output_content = starfit(input_contents, z_max=30,
            z_exclude=("Li", ), combine_CNO=combine_CNO)

        with open("{0}/{1}.txt".format(RESULTS_FOLDER, prefix), "w") as fp:
            fp.write(input_contents)
        with open("{0}/{1}.pkl".format(RESULTS_FOLDER, prefix), "wb") as fp:
            pickle.dump((result, payload, output_content), fp, -1)
        with open("{0}/{1}.html".format(RESULTS_FOLDER, prefix), "w") as fp:
            fp.write(output_content)


    # Perturb stars within the uncertainties
    for star in stars:
        for i in range(N_PERTURBATIONS):

            print("Star {0} perturbation {1}".format(star, i))

            prefix = clean_star_name(star)

            prefix_filename = "{0}/{1}.perturb-{2}".format(RESULTS_FOLDER, prefix, i)
            if os.path.exists("{0}.pkl".format(prefix_filename)) and not OVERWRITE:
                print("Skipping {0} perturbation {1} because {2}.pkl exists..".format(star, prefix, prefix_filename))
                continue

            try:
                input_contents, elements, upper_limits, z_max, combine_CNO = starfit_format(star, perturb=True)

            except AssertionError:
                print("Skipping {}..".format(star))
                continue

            except:
                raise

            result, payload, output_content = starfit(input_contents,
                z_max=30, z_exclude=("Li", ), combine_CNO=combine_CNO)

            with open("{0}.txt".format(prefix_filename), "w") as fp:
                fp.write(input_contents)
            with open("{0}.pkl".format(prefix_filename), "wb") as fp:
                pickle.dump((result, payload, output_content), fp, -1)
            with open("{0}.html".format(prefix_filename), "w") as fp:
                fp.write(output_content)

    raise a
    """

    raise a

    # Get RA, DEC from stars:



    raise a


    #result, payload, content = starfit(os.path.abspath("data/example_keller_star.dat"))

    """


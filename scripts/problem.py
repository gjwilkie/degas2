import species 
#import PMI as pm
#import material as mat
#import reaction as rc
import utils
import os
import netCDF4 as nc
import numpy as np

class Problem:
    """
    A python representation of a DEGAS2 'problem'. Contains all the data needed
    to generate a problem_infile.

    Attributes:
        description: Arbitrary string representing the problem name. (optional)
        test_species: List of strings specifying test particle species. Starts with '0' for 'ghost' geometry species
        background_species: List of strings specifying background plasma species. By convention, first element is 'e' for electrons.
        reactions: List of strings specifying each reaction in the problem. See data/reactions.input
        materials: List of strings specifying each material handled in the problem. See data/materials.input
        pmi: List of strings specifying each PMI handled in the problem. See data/pmi.input
    """
    def __init__(self):
        """
        Initializes a problem. Takes no arguments.
        """
        self.description = "No description"
        self.test_species = []
        self.background_species = []
        self.reactions = []
        self.materials = []
        self.pmi = []

    def add_test_species(self,sp):
        """
        Appends a string to the problem's test_species attribute

        Args:
            test: Species string
        """
        self.test_species.append(sp)

    def add_background_species(self,sp):
        """
        Appends a string to the problem's background_species attribute

        Args:
            sp: Species string
        """
        self.background_species.append(sp)

    def add_reaction(self,reaction):
        """
        Appends a string to the problem's reactions attribute

        Args:
            reaction: String for the reaction label in data/reactions.input
        """
        self.reactions.append(reaction)

    def add_PMI(self,pmi):
        """
        Appends a string to the problem's pmi attribute

        Args:
            pmi: String for the PMI label in data/pmi.input
        """
        self.pmi.append(pmi)

    def add_material(self,material):
        """
        Appends a string to the problem's materials attribute

        Args:
            material: String for the PMI label in data/materials.input
        """
        self.materials.append(material)

    def parse_input_file(self,filename):
        """
        Reads a degas2 problem_infile and stores its contents to class instance self.

        Args:
            filename: String for the relative path to the file from which to read in problem.
        """
        headers = ["TEST","BACKGROUND","REACTION","MATERIALS","PMI"]
        d = utils.parse_file_with_headers(filename,headers)
        for name in d["TEST"]:
            self.add_test_species( species.lookup(name) )
        for name in d["BACKGROUND"]:
            self.add_background_species( species.lookup(name) )
        for name in d["REACTION"]:
            self.add_reaction(reaction.lookup(name) )
        for name in d["MATERIALS"]:
            self.add_material(material.lookup(name) )
        for name in d["PMI"]:
            self.add_PMI(pmi.lookup(name) )

    def generate_input_file(self,filename):
        """
        Generates a degas2 problem_infile from class instance

        Args:
            filename: String for the relative path to the file to write out the problem data
        """
        d={}
        d["TEST"] = self.test_species
        d["BACKGROUND"] = self.background_species
        d["REACTION"] = self.reactions
        d["MATERIALS"] = self.materials
        d["PMI"] = self.pmi
        utils.write_file_with_headers(filename,d,"$ $")
   
def defineProblem(testSps,backSps,reactions,mats,pmis):
    """
    Generates a degas2 problem with one command.

    Args:
        testSps: Array of strings for test species. First element must be '0'
        backSps: Array of strings for background species. First element is 'e' by convention if present.
        reactions: Array of strings for reactions to include.
        mats: Array of strings for materials to include.
        pmis: Array of strings for PMI to include.

    Returns:
        A Problem with the specified properities.
    """

    p = Problem()
    for t in testSps:
        p.add_test_species(t)
    for b in backSps:
        p.add_background_species(b)
    for r in reactions:
        p.add_reaction(r)
    for m in mats:
        p.add_material(m)
    for m in pmis:
        p.add_PMI(m)
    return p

def generateProblemInput(testSps,backSps,reactions,mats,pmis,filename="problem.in"):
    """
    Generates a degas2 problem with one command and writes the input file

    Args:
        testSps: Array of strings for test species. First element must be '0'
        backSps: Array of strings for background species. First element is 'e' by convention if present.
        reactions: Array of strings for reactions to include.
        mats: Array of strings for materials to include.
        pmis: Array of strings for PMI to include.
        filename: String for the relative path of the input filename (optional)

    Returns:
        A Problem with the specified properities.
    """
    p = defineProblem(testSps,backSps,reactions,mats,pmis)
    p.generate_input_file(filename)
    return p

def genStdProblem(label):
    """
    Generates a degas2 problem input using a single string for a select 'standard' configurations.

    Args:
        label: String. One of: "C-H" (hydrogen against carbon walls - no molecules), "C" (hydrogen against carbon walls including molecules), "Li" (hydrogen against lithium walls), "Li_reflOnly" (hydrogen on lithium with no desorptioon), "LiOH" (lithium against LiOH including molecules), "Fe_and_mirror" (iron and mirror including molecules).

    Returns:
        A Problem of the specified standard type.
    """
    if label == "C-H":
        p=generateProblemInput(["0","H"],["e","H+"],["hionize5","hh_chargex"],["C"],["hdesorbc_xgc","hreflc"])
    if label == "C-D-xgc":
        p=generateProblemInput(["0","D","D2","D2+"],["e","D+"],["hionize5","dchex_const","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec","hrecombine5"],["C"],["hdesorbc","h2desorbc","dreflc"])
    if label == "C":
        p=generateProblemInput(["0","H","H2","H2+"],["e","H+"],["hionize5","hh_chargex","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec"],["C"],["hdesorbc","h2desorbc","hreflc"])
    if label == "C-D":
        p=generateProblemInput(["0","D","D2","D2+"],["e","D+"],["hionize5","dd_chargex","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec"],["C"],["hdesorbc","h2desorbc","dreflc"])
    if label == "C-D2":
        p=generateProblemInput(["0","D","D2","D2+"],["e","D+"],["hionize5","hrecombine5","dd_chargex","h2dis_jcr","h2ion_jcr","h2dision_jcr","h2pdision_jcr","h2pdis_jcr","h2pdisrec_jcr"],["C"],["hdesorbc","h2desorbc","dreflc"])
    if label == "mirror-D":
        p=generateProblemInput(["0","D","D2","D2+"],["e","D+"],["hionize5","dd_chargex","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec"],["mirror"],["hmirror","h2mirror"])
    if label == "Li":
        p=generateProblemInput(["0","H"],["e","H+"],["hionize5","hh_chargex"],["Li"],["H_refl_svftrim_Li","hdesorbLi"])
    if label == "Li_reflOnly":
        p=generateProblemInput(["0","H"],["e","H+"],["hionize5","hh_chargex"],["Li"],["H_refl_svftrim_Li"])
    if label == "LiOH":
        p=generateProblemInput(["0","H","H2","H2+"],["e","H+"],["hionize5","hh_chargex",\
                "h2dis","h2dis_n3","h2ion","h2dision","h2pdision","h2pdis"],["LiOH"],\
                ["H_refl_vftrim_LiOH","h_des_maxw_LiOH","h2_des_maxw_LiOH"])
    if label == "Fe_and_mirror":
        p=generateProblemInput(["0","H","H2","H2+"],["e","H+"],["hionize5","hh_chargex","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec"],["Fe","mirror"],["h_des_maxw_fe","h2_des_maxw_fe","hreflfe","hmirror","h2mirror"])
    return p

def get_reactions_from_problem(degas2_in, problem_nc, verbose=True):
    """ Function used to return information about reactions in a given degas2 simulation.
    This is designed to read,
      (1) the degas2.in file to figure out where the reactions.nc file is...
      (2) the problem.nc file to figure out what reactions were requested by the user. 

    Args:
        degas2_in (str) path to the degas2.in file.
        problem_nc (str) path to the problem_nc file.

    Returns:
        dictionary with a subset of the info from reactions.input pertaining to this problem.nc.
    """
    for f in [degas2_in, problem_nc]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"The file {f} does not exist.")
    # Read degas2_in to determine the reactions.input file path,
    with open(degas2_in, 'r') as file:
        dg2in_lines = file.readlines()
        #reactions_in = [l.strip() for l in dg2in_lines if l.strip().startswith("reaction_infile")][0].split(" ")[1]
        reactions_nc = [l.strip() for l in dg2in_lines if l.strip().startswith("reactionfile")][0].split(" ")[1]
    
    # use helper function to prase the reactions.input file,
    #reactions_dict = read_reactions_input(reactions_in)
    reactions_dict = read_reactions_nc(reactions_nc)

    # Read the problem.nc file and figure out the indices of reactions used in this problem.
    # This *appears* to be held in the 'problem_rc' data field.
    with nc.Dataset(problem_nc, 'r') as ds:
        problem_rc = ds.variables['problem_rc'][:] # np.MaskedArray
        rxn_in_problem = [ list(reactions_dict.keys())[i-1] for i in problem_rc ]
    out = {k: reactions_dict[k] for k in rxn_in_problem}
    return out

def read_reactions_nc(reactions_nc):
    """ Function to prase the reactions.nc data file.
    This function is similar to problem.read_reactions_input() but it reads directly from 
    the NetCDF file rather than the .input text file.
    
    Args:
        reactions_nc (str) path to the reactions.nc file.

    Returns:
        dictionary with info about reactions included in DEGAS2.   
    """
    if not os.path.exists(reactions_nc):
        raise FileNotFoundError(f"The file {reactions_nc} does not exist.")

    # Open the NetCDF file in read mode
    with nc.Dataset(reactions_nc, 'r') as ds:
        # helper to safely decode byte arrays to strings
        def decode_strings(var_name):
            data = ds.variables[var_name][:]
            if data.ndim > 1:
                data = nc.chartostring(data)
            # NetCDF4 usually returns numpy object arrays for strings, 
            # but we ensure they are standard strings for the dict keys/values
            out = [str(s).rstrip() if not isinstance(s, bytes) else s.decode('utf-8') for s in data]
            return out

        # Extracting the variables
        data_path = decode_strings('reaction_filename')
        long_name = decode_strings('reaction_name')
        classification = decode_strings('reaction_type')
        name = decode_strings('reaction_sy')

    # Constructing the dictionary-of-dictionaries
    result = {}
    for i in range(len(name)):
        k = name[i]
        result[k] = {
            'data_path': data_path[i],
            'long_name': long_name[i],
            'classification': classification[i]
        }
        
    return result
    
def read_reactions_input(reactions_input):
    """ Function used to read the reactions.input file and convert it to a dict-of-dicts
    Ignores lines starting with '#' or '$'.
    """
    result = {}
    
    if not os.path.exists(reactions_input):
        raise FileNotFoundError(f"The file {reactions_input} does not exist.")

    with open(reactions_input, 'r') as f:
        # 1. Filter out comments and empty lines immediately
        lines = [
            line.strip() for line in f 
            if line.strip() and not line.startswith(('#', '$'))
        ]

    # 2. Process the cleaned lines in steps of 4
    # Format: 0:long_name, 1:name+symbol, 2:data_path, 3:classification
    for i in range(0, len(lines), 4):
        try:
            full_name = lines[i]
            
            # Split 'hionize e + H * -> e + H+ + e' into ['hionize', 'rest...']
            # maxsplit=1 ensures we only split at the first whitespace/tab
            name_parts = lines[i+1].split(None, 1)
            name = name_parts[0]
            symbolic_rep = name_parts[1] if len(name_parts) > 1 else ""
            
            data_path = lines[i+2]
            classification = lines[i+3]

            # 3. Build the nested dictionary
            result[name] = {
                "long_name": full_name,
                "symbolic_representation": symbolic_rep,
                "data_path": data_path,
                "classification": classification
            }
        except IndexError:
            print(f"WARN: Found a partial text block at the end of the file starting with '{lines[i]}'. Skipping.")
            break

    return result

if __name__ == "__main__":
    # Run a test of some functions in this file,
    d2path = os.environ["DEGAS2_DIR"] # set with >> module load degas2
    print(f"TEST (problem.py): DEGAS2_DIR={d2path}")

    rxn_nc_path = "/data/reactions.nc"
    test = read_reactions_nc(d2path + rxn_nc_path)
    print(f"TEST (problem.py): Found {len(test.keys())} reactions in $DEGAS2_DIR{rxn_nc_path}")
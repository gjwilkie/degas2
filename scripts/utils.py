import netCDF4 as nc

# TODO
# Returns a dictionary with keys as "mandatory_headers"
def parse_file_with_headers(filename,mandatory_headers):
    f = open(filename,"r")
    lines = f.readlines()
    f.close()

    #def get_header_indices(lines,headers):

    idx = get_header_indices(lines,mandatory_headers)

    d = {}

    utils.write_file_with_headers(filename,d,"$ $")

def write_file_with_headers(filename,d,info=None):
    """
    Writes a simple degas2 input file using an input dictionary.
    Currently only works for problemsetup input files

    Args:
        filename: the name of the file to write
        d: dictionary of data where the key is a header for the group. Individual objects within each key are strings which are written line-by-line.
    """

    f = open(filename,"w")
    if info:
        f.write(info+"\n")
    for key in d:
        f.write(key+"\n")
        for obj in d[key]:
            f.write(obj+"\n")
    f.close()


def get_strings_1d(self,s):
    data = []
    for i in range(0,len(self[s])):
         data.append(str(nc.chartostring(self[s])).strip())
    return data


def str_to_intarr(s, N=None):
    """ Helper to convert special strings to integer arrays.
    This is used in the source.py:Source class to help define source segments/strata.
    Examples: 
        '*'       --> [0, 1, ... N-1] (if N=None --> '*')
        '1 2 10'  --> [1, 2, 10]
        '1:3'     --> [1, 2]
        '1:3+7:9' --> [1, 2, 7, 8]
    :arg s: string to parse
    :arg N: number of segments in total - if None, may return '*'.
    """
    # 1. Input validation for allowed characters
    # Note: Spaces are explicitly allowed to support the '1 2 10' format
    allowed_chars = set("0123456789*: +")
    if not set(s).issubset(allowed_chars):
        raise ValueError("Invalid characters in string. Only numbers, '*', ':', '+', and spaces are permitted.")
    
    # 2. Parse the string
    # Replace '+' with spaces to treat all separated groups uniformly
    normalized_str = s.replace('+', ' ')
    tokens = normalized_str.split()
    
    result = []
    
    for token in tokens:
        if token == '*':
            if N is None:
                 # if the user doesn't tell us the total number, we don't know how to handle the wildcard.
                return '*'
            else:
                result.extend(range(N))
        elif ':' in token:
            parts = token.split(':')
            if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isdigit():
                raise ValueError(f"Invalid range format: '{token}'")
            start, end = int(parts[0]), int(parts[1])
            result.extend(range(start, end))
        else:
            if not token.isdigit():
                raise ValueError(f"Invalid segment format: '{token}'")
            result.append(int(token))
    
    if N is not None:
        # 3. Validate against total segment count (N)
        # Assuming standard 0-indexed arrays, values >= N are out of bounds
        if any(val >= N for val in result):
            raise ValueError(f"Segment number cannot be greater than or equal to N ({N}).")
        
    # Optional: Deduplicate values and return as a sorted list
    return sorted(list(set(result)))
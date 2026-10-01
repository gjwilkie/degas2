from polygon import *
from utils import str_to_intarr
class Source:
    """ Class managing DEGAS2 sources.
    """
    def __init__(self,nflights,stype,species,rootspecies=None,specify_flux=True,sourcefile="sourcefile.txt",
                 sourcefile_fmt="row",pufftemp=None,puffexp=None,strength=None,stratum=None,segment=None,
                 e_bin_num=10,e_bin_min=0.1,e_bin_max=1000.0,e_bin_log=True):
        """ Initialization of DEGAS Source class.
            Args:
                nflights: (int) Number of flights
                stype: (str) one of ["plate", "puff", "vol_source", "plt_e_bins", "snapshot", "recomb"]
                species: (str) species /generated/ by the source, must be in list of test species (e.g. D)
            Kwargs:
                rootspecies: (str) species underlying the source (e.g. D+)
                specify_flux: (bool) informs the code that 'source_strength' is a flux rather than current.
                sourcefile: (str) filename for tabular file-based input.
                sourcefile_fmt: (str) one of ['tabular','row']
                pufftemp: (float) Temperature of the puff in Kelvin.
                strength: (list of float) strength of the source, units are determined by specify_flux.
                           The defineback program expects one source_strength for every stratum, segment.
                           The values are separated by spaces: S1 S2 S3...
                stratum: (list of int) Used with strength and segment to specify the location.
                                The values are separated by spaces: 3 3 3
                segmment: (list of int) Used with strength, stratum to specify the location of the source.
                                The values are separated by spaces: 100 200 300
        """
        self.nflights = nflights
        self.type = stype
        self.species = species
        # Sanitize input,
        if stype == "vol_source":
            strength=None
            sourcefile_fmt='row'
            self.geom = 'volume'
            pufftemp=None
            puffexp=None
        else:
            self.geom = 'surface'
        #
        self.sourcefile_fmt=sourcefile_fmt
        #
        if rootspecies is None:
            # automatically configure based on stype...
            if stype == "plate":
                rootspecies = species + "+"
            else:
                rootspecies = species
            self.rootspecies = rootspecies
        else:
            self.rootspecies = rootspecies
        #
        if specify_flux:
            self.specify_units = "specify_flux"
        else:
            self.specify_units = "specify_current"
        #
        if pufftemp:
            self.pufftemp = pufftemp
        else:
            self.pufftemp=None
        #
        if puffexp:
            self.puffexp = puffexp
        else:
            self.puffexp = None
        # 
        if strength:
            self.sourcefile = None
            # The source is to be specified with a combination of 'strength', 'stratum', and 'segment'
            # However, if the given strength/segment/stratum lists are too long, this method will write a sourcefile.
            self.set_source_from_segments(strength, stratum, segment)
            self.check_segments(sourcefile=sourcefile)
        else:
            self.sourcefile = sourcefile

        if stype=="plt_e_bins":
            self.type = "plt_e_bins"
            self.e_bin_num = e_bin_num
            self.e_bin_min = e_bin_min
            self.e_bin_max = e_bin_max
            if e_bin_log:
                self.e_bin_spacing = "log"
            else:
                self.e_bin_spacing = "linear"
            
    def set_source_from_segments(self, strength, stratum, segment):
        """ Method to set Source attrs given: strength, stratum, segment.
        Args:
            strength: (str, float, or list-of-float)
                      The correct attr is created if strength=1e24, '1e24' or a list thereof
            stratum: (str, int, or list-of-int)
                     Can be a special formatted string parsed by `str_to_intarr` (e.g. "10:20").
                     If it resolves to a single int, it will be distributed across segments.
                     Cannot be the '*' wildcard.
            segment: (str, int, or list-of-int)
                     Can be a special string (e.g. "100:200", "100 200", "100:200+10:20").
                     "*" = all segments wildcard, strength and stratum evaluate to a single value.
        """
        # 0. Strict validation for stratum
        if stratum == "*" or (isinstance(stratum, list) and "*" in stratum):
            raise ValueError("Stratum cannot be the '*' wildcard.")

        # 1. Parse segment string (or handle early if it's the wildcard)
        if isinstance(segment, str):
            segment_list = str_to_intarr(segment, N=None) # handles special strings.
        else:
            segment_list = segment if isinstance(segment, list) else [segment]
            segment_list = [int(s) for s in segment_list] # ensure user has provided ints.

        # 2. Handle wildcard "*" globally
        # This triggers either if segment was passed as "*" or if str_to_intarr fell through and returned "*"
        if segment_list == "*":
            if not isinstance(strength, list): strength = [strength]
            if not isinstance(stratum, list): stratum = [stratum]
            
            self.segment = "*"
            self.strength = f"{float(strength[0]):.3e}"
            self.stratum = str(int(stratum[0]))
            return

        # 3. Parse stratum into lists of ints
        if isinstance(stratum, str):
            stratum_list = str_to_intarr(stratum, N=None) # handles special strings.
        else:
            stratum_list = stratum if isinstance(stratum, list) else [stratum]
            stratum_list = [int(s) for s in stratum_list]

        # 4. Distribute/broadcast arrays to match lengths
        len_seg = len(segment_list)
        len_strat = len(stratum_list)

        if len_seg == 1 and len_strat > 1:
            segment_list = segment_list * len_strat
        elif len_strat == 1 and len_seg > 1:
            stratum_list = stratum_list * len_seg
        elif len_seg != len_strat:
            raise ValueError(f"Shape mismatch: segment has {len_seg} elements, stratum has {len_strat}.")

        N = len(segment_list)

        # 5. Handle strength normalization (Flux vs Current)
        if not isinstance(strength, list):
            strength = [strength]
            
        strength_list = [float(s) for s in strength]
        if len(strength_list) == 1:
            # Default to specify_current behavior (division) if specify_units is missing/different
            units = getattr(self, "specify_units", "specify_current")
            if units == "specify_flux":
                # Repeat the given flux value over the segments
                strength_list = [strength_list[0]] * N
            else:
                # Divide the given total 'strength' over the number of segments
                strength_list = [strength_list[0] / N] * N
        elif len(strength_list) == N:
            pass
        else:
            raise ValueError(f"Shape mismatch: strength has {len(strength_list)} elements, segments/strata has {N}.")

        # 6. Final assignment to string attributes expected by DEGAS2
        self.segment = " ".join(f"{s}" for s in segment_list)
        self.stratum = " ".join(f"{s}" for s in stratum_list)
        self.strength = " ".join(f"{s:.6e}" for s in strength_list)

    def check_segments(self, max_line_len=275, sourcefile="sourcefile.txt"):
        """ Method to check the segment/strength/stratum data and reformat if needed into a sourcefile.
        Lines in the 'db.in' file (below) can only be LINELEN(=300 for this build), if they're too long the
        Source must be specified via a sourcefile.txt instead of directly in db.in.
        """
        line_lens = [len(self.segment), len(self.strength), len(self.stratum)]
        if any([l >= max_line_len for l in line_lens]):
            print(f"WARNING: Some lines in 'db.in' may exceed the maximum line length! - writing equivalent '{sourcefile}'")
            # Convert the segment, strength, and stratum back into lists...
            strength = [float(s) for s in self.strength.split(" ")]
            stratum = [int(s) for s in self.stratum.split(" ")]
            segment = [int(s) for s in self.segment.split(" ")]
            write_sourcefile_from_segments(strength=strength, stratum=stratum, segment=segment, filename=sourcefile)
            self.sourcefile = sourcefile # This will use the sourcefile instead of the strength/stratum/segment attrs. 
            self.sourcefile_fmt = "tabular" 
            
        return

def write_sourcefile_from_segments(strength=[1e24], stratum=[3], segment=[1], filename="sourcefile.txt"):
    """ Function to write a sourcefile.txt in TABULAR format defining more distributed sources.
    Kwargs:
        strength:
        stratum:
        segment:
        filename:
    """
    f = open(filename, "w")
    f.write("# Sourcefile automatically generated by script\n")
    f.write("# i_stratum i_segment S_strength\n")
    for strat, seg, S in zip(stratum, segment, strength):
        f.write(f"{int(strat)} {int(seg)} {S}\n")
    f.close()
    
def write_db_input(source_groups,plasmafile="plasmafile.txt",filename="db.in",init=False,t0=-1.0,tf=-1.0):
    """ Function to write the 'db.in' input file for the defineback program.
    String formatting is important here because these values will be read by FORTRAN
    Args:
        source_groups: (Source type) object with necessary data as attrs.
        plasmafile: (str) filename containing the plasma background.
        filename: (str) filename written by this function.

    WARNING: in defineback.web, lines from this file are read into the variable 'line'
        declared with,
            >> character*LINELEN, line (e.g. ln 1088 in defineback.web)
        Thus, there is an implicit maximum line length limiting the number of sources.
        This common variable "LINELEN" appears to be set to 300 in 'macros.hweb'.
    """
    Nsource=len(source_groups)
    f = open(filename,"w")
    if tf > 0:
        if init:
            f.write("time_initialization\n")
        f.write("time_interval %8.3e %8.3e\n"%(t0,tf)) 
    f.write(f"plasma_file {plasmafile}\n")
    # For each source group, 
    for i in range(0, Nsource):
        source = source_groups[i]
        f.write("new_source_group\n")
        f.write(f"  source_type {source.type}\n")
        f.write(f"  source_geom {source.geom}\n")
        f.write(f"  source_species {source.species}\n")
        f.write(f"  source_root_sp {source.rootspecies}\n")
        f.write(f"  {source.specify_units}\n")
        # if sourcefile exists, stratum/segment/strength info is included there,
        if source.sourcefile is not None:
            f.write(f"  source_file {source.sourcefile} {source.sourcefile_fmt}\n")
        else:
            f.write("  source_stratum "+str(source.stratum)+"\n")
            f.write("  source_segment "+str(source.segment)+"\n")
            f.write("  source_strength "+str(source.strength)+"\n")

        if source.type == "puff":
            if source_groups[i].pufftemp != None:
                f.write("  source_puff_temp "+str(source.pufftemp)+"\n")
            if source_groups[i].puffexp != None:
                f.write("  source_puff_exponent "+str(source.puffexp)+"\n")
        
        if source_groups[i].type == "plt_e_bins":
            f.write("  source_e_bin_min "+str(source.e_bin_min)+"\n")
            f.write("  source_e_bin_max "+str(source.e_bin_max)+"\n")
            f.write("  source_e_bin_num "+str(source.e_bin_num)+"\n")
            f.write("  source_e_bin_spacing "+source.e_bin_spacing+"\n")

        f.write("  source_nflights "+str(source.nflights)+"\n")
        f.write("end_source_group\n")
    f.close()



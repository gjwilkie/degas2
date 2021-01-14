# This contains all the information needed for the datasetup executable. These nc files live in the data directory, and rarely change from run to run, unless new data is added.
import netCDF4 as nc
import sys

class degas2_data:
    def __init__(self,datadir):
        self.dir = datadir

        self.elements_infile, self.elementsfile = self.getfiles("elements")
        self.species_infile, self.speciesfile = self.getfiles("species")
        self.reaction_infile, self.reactionfile = self.getfiles("reactions")
        self.materials_infile, self.materialsfile = self.getfiles("materials")
        self.pmi_infile, self.pmifile = self.getfiles("pmi")

    def getfiles(self,filetype):
        infile = self.dir + "/" + filetype + ".input"
        outfile = self.dir + "/" + filetype + ".nc"
        return infile, outfile

    #def writefile(self):

    #def run(self):

# Some utility functions:

# Accepts a raw netcdf variable representing an array of strings.
# Returns an array of stripped strings as native python strings.
def convertStr(ncvar):
    ndims=len(ncvar.shape)
    if ndims == 1:
        return str(nc.chartostring(ncvar)).strip()
    elif ndims ==2:
        result = []
        for i in range(0,ncvar.shape[0]):
            result.append( str(nc.chartostring(ncvar[i])).strip() )
        return result
    elif ndims ==3:
        result = []
        for i in range(0,ncvar.shape[0]):
            temp1 = []
            for j in range(0,ncvar.shape[1]):
                temp1.append(str(nc.chartostring(ncvar[i,j])).strip())
            result.append(temp1)
        return result
    elif ndims ==4:
        result = []
        for i in range(0,ncvar.shape[0]):
            temp1 = []
            for j in range(0,ncvar.shape[1]):
                temp2 = []
                for k in range(0,ncvar.shape[2]):
                    temp2.append(str(nc.chartostring(ncvar[i,j,k])).strip())
                temp1.append(temp2)
            result.append(temp1)
        return result
    else:
        sys.exit("Can only handle up to 3-dimensional string arrays. Add a new case for dimensionality "+str(ndims))
        
# Class containing interpolant function to get atomic physics data from a degas2 netcdf data file.
#   filename: a string for the path and file name of the netcdf file to query
#   requestedData: a string representing the data to read (usually "cross_section" or "reaction_rate")
class d2_atomic_table:
    def __init__(self,filename,requestedData):
        data=nc.Dataset(filename)
        self.dataName = requestedData
        self.description = data.getncattr("data_version")

        varNames = convertStr(data["xs_var"])
        table_idx = varNames[:][0].index(requestedData) + 1
        self.rank = data["xs_rank"][table_idx]
        self.indepVarNames = varNames[table_idx][1:]
        self.dims = data["xs_tab_index"][table_idx][0:self.rank]

        xs_spacing = convertStr(data["xs_spacing"][table_idx][0:self.rank+1])
        
        if xs_spacing[0] == "log":
            self.indepLog = True

        self.depLog = [False]*self.rank
        for i in range(0,self.rank):
            if xs_spacing[i+1] == "log":
                self.indepLog[i] = True

        xs_units = convertStr(data["xs_units"][table_idx][0:self.rank+1])
        self.depUnits = xs_units[0]
        self.indepUnits = xs_units[1:self.rank+1]
        
        self.depMult = data["xs_mult"][table_idx][0]
        self.indepMult = data["xs_mult"][table_idx][1:self.rank+1]

        inddepTable = np.zeros(self.dims)

        self.indepMin = self.indepMult[:] * data["xs_min"][table_idx][0:self.rank]
        self.indepMax = self.indepMult[:] * data["xs_min"][table_idx][0:self.rank]

        self.delta = np.zeros(self.rank)
        for i in range(0,self.rank):
            if indepLog[i] == True:
                self.indepMin[i] = np.log(indepMin[i])
                self.indepMax[i] = np.log(indepMax[i])
            self.delta[i] = (indepMax[i]-indepMin[i])/(self.dims[i]-1)

        base_idx = data["xs_data_base"][table_idx] + 1
        ndata = data["xs_data_inc"][table_idx] 
        self.depTable = self.depMult * data["xs_data_tab"][base_idx:base_idx + ndata]
        depTable.reshape(self.dims)
        if depLog == True:
            indepTable = np.log(indepTable)

    # Queries a table already created with an array for independent variables. Linearly interpolates on table.
    def getData(self,indepVals_in):

        indepVals = np.copy(indepVals_in)

        # Get the lower index for interpolation
        idx_lower = []
        idx_upper = []
        weights = []
        for i in range(0,self.rank):
            if self.indepLog[i]:
                indepVals[i] = np.log(indepVals_in[i])

            idx.append( min( self.dims[i], \
                    int( (indepVals[i] - self.indepMin[i]) //self.delta[i] ) ) )

            weights.append( [ \
                self.indepMin[i] + self.delta[i]*(idx[i]+1) - indepVals[i]  , \
                indepVals[i] - (self.indepMin[i] + self.delta[i]*(idx[i]) ) ] )
            weights[i][:] = weights[i][:]/self.delta[i]

        if self.rank == 0:
            result = self.depTable[0]
        elif self.rank == 1:
            result = weights[0][0]*(self.indepMin[0] + self.delta[0]*idx[i]) + \
                    weights[0][1]*(self.indepMin[0] + self.delta[0]*(idx[i]+1)) 
        elif self.rank == 2:
            result = 0.0
            for i in range(0,2):
                for j in range(0,2):
                    result += weights[i][0]*weights[j][1]*self.depTable[idx[i],idx[j]]
        elif self.rank == 3:
            result = 0.0
            for i in range(0,2):
                for j in range(0,2):
                    for k in range(0,2):
                        result += weights[i][0]*weights[j][1]*weights[k][2]*self.depTable[idx[i],idx[j],idx[k]]
        else:
            sys.exit("Routine only handles up to rank 3 interpolation. Update the code.")

        if depLog:
            result = np.exp(result)

        return result

    # Prints to stdout various information about the table
    def printInfo(self):
        print("Dependent variable: "+self.dataName+" for process "+self.description+"\n")
        print("All input and output units are SI.\n")
        print("Extrapolations beyond table bounds are nearest-neighbor.\n")
        if self.indepLog):
            print("Logarithmic scale\n")
        else:
            print("Linear scale\n")
        print("Rank = %d \n" % self.rank)
        print("Independent variables:\n")
        for idim in range(0,self.rank):
            print("  "+self.indepVarName[idim]+":\n")
            print("    Number of datapoints: %d \n" % self.dims[i] )
            print("    Table minimum: %f "+ indep_var_units_[i]+"\n" % (self.indepMin[i]/self.indepMult[i]) )
            print("    Table maximum: %f "+ indep_var_units_[i]+"\n" % (self.indepMax[i]/self.indepMult[i]) )
            if self.indepLog[i]:
                print("    Logarithmic scale\n")
            else:
                print("    Linear scale\n")




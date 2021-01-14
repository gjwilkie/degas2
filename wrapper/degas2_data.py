# This contains all the information needed for the datasetup executable. These nc files live in the data directory, and rarely change from run to run, unless new data is added.
import netCDF4 as nc
import sys
import numpy as np
from importlib import reload

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
        depNames = [varNames[i][0] for i in range(0,len(varNames))]
        
        table_idx = depNames.index(requestedData)
        self.rank = int(data["xs_rank"][table_idx])
        self.indepVarNames = varNames[table_idx][1:]
        self.dims = np.array(data["xs_tab_index"][table_idx][0:self.rank])

        xs_spacing = convertStr(data["xs_spacing"][table_idx][0:self.rank+1])
        
        if xs_spacing[0] == "log":
            self.depLog = True

        self.indepLog = [False]*self.rank
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
        self.indepMax = self.indepMult[:] * data["xs_max"][table_idx][0:self.rank]

        self.delta = np.zeros(self.rank)
        for i in range(0,self.rank):
            if self.indepLog[i] == True:
                self.indepMin[i] = np.log(self.indepMin[i])
                self.indepMax[i] = np.log(self.indepMax[i])
            self.delta[i] = (self.indepMax[i]-self.indepMin[i])/(self.dims[i]-1)

        base_idx = data["xs_data_base"][table_idx] 
        ndata = data["xs_data_inc"][table_idx] 
#        self.depTable = np.array(self.depMult * data["xs_data_tab"][base_idx:(base_idx + ndata)])#.reshape(self.dims[-1::-1])
        self.depTable = np.array(self.depMult * data["xs_data_tab"][base_idx:(base_idx + ndata)]).reshape(self.dims[-1::-1])
        self.depTable = np.transpose(self.depTable)

        if self.depLog == True:
            self.depTable = np.log(self.depTable)

    # Queries a table already created with an array for independent variables. Linearly interpolates on table.
    def getData(self,indepVals_in):

        indepVals = np.array(indepVals_in)

        # Get the lower index for interpolation
        idx = []
        weights = []
        for i in range(0,self.rank):
            if self.indepLog[i]:
                indepVals[i] = np.log(indepVals[i])

            if indepVals[i] < self.indepMin[i]:
                idx.append(0)
                weights.append([1.0,0.0])
            elif indepVals[i] > self.indepMax[i]:
                idx.append(self.dims[i]-2)
                weights.append([0.0,1.0])
            else:
                idx.append( min( self.dims[i], \
                        int( (indepVals[i] - self.indepMin[i]) //self.delta[i] ) ) )

                weights.append( [ \
                    self.indepMin[i] + self.delta[i]*(idx[i]+1) - indepVals[i]  , \
                    indepVals[i] - (self.indepMin[i] + self.delta[i]*(idx[i]) ) ] )
                weights[i][:] = weights[i][:]/self.delta[i]

        print(idx)
        print(weights)
        if self.rank == 0:
            result = self.depTable[0]
        elif self.rank == 1:
            result = weights[0][0]*self.depTable[idx[0]] + \
                    weights[0][1]*self.depTable[idx[0]+1]
        elif self.rank == 2:
            result = 0.0
            for i in range(0,2):
                for j in range(0,2):
                    result += weights[0][i]*weights[1][j]*self.depTable[idx[0]+i,idx[1]+j]
        elif self.rank == 3:
            result = 0.0
            for i in range(0,2):
                for j in range(0,2):
                    for k in range(0,2):
                        result += weights[0][i]*weights[1][j]*weights[2][k]*self.depTable[idx[0]+i,idx[1]+j,idx[2]+k]
        else:
            sys.exit("Routine only handles up to rank 3 interpolation. Update the code.")


        if self.depLog:
            result = np.exp(result)

        return result

    # Prints to stdout various information about the table
    def printInfo(self):
        print("Dependent variable: "+self.dataName+" for process "+self.description)
        print("All input and output units are SI.")
        print("Extrapolations beyond table bounds are nearest-neighbor.")
        if self.indepLog:
            print("Logarithmic scale")
        else:
            print("Linear scale")
        print("Rank = %d " % self.rank)
        print("Independent variables:")
        for i in range(0,self.rank):
            print("  "+self.indepVarNames[i]+":")
            if self.indepLog[i]:
                minval = np.exp(self.indepMin[i])/self.indepMult[i]
                maxval = np.exp(self.indepMax[i])/self.indepMult[i]
            else:
                minval = self.indepMin[i]/self.indepMult[i]
                maxval = self.indepMax[i]/self.indepMult[i]
            print("    Number of datapoints: %d " % self.dims[i] )
            print("    Table minimum: %e %s" % (minval,  self.indepUnits[i]) )
            print("    Table maximum: %e %s" % (maxval,  self.indepUnits[i]) )
            if self.indepLog[i]:
                print("    Logarithmic scale")
            else:
                print("    Linear scale")




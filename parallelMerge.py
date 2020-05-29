import json, sys, gzip, time, os
import multiprocessing as mp
from multiprocessing import Pool, cpu_count


class DataAnalyzer:
    """Parallel processing data and generate mRNA Matrix."""
    def __init__(self, path, outFileName):
        # metadata
        with open(path, 'r') as jf:
            self._data = json.load(jf)
            
        # use to count program calculation time
        self._startTime = time.time()
        # final result
        self._results = {}
        self._normalSamples = []
        self._tumorSamples = []
        self._outFileName = outFileName
        
    def run(self):
        """Separate task depends on number of cores"""
        n = cpu_count() if len(self._data) >= cpu_count() else 1
        
        pool = Pool(processes=n)
        poolResults = []
        j = 0
        k = len(self._data) // n
        # separate task then distribute to each process
        for i in range(n):
            if i == n - 1:
                data = self._data[j:]
            else:
                data = self._data[j:j+k]
                j += k
            poolResults.append(pool.map_async(self._getResult, (data,)))
            
        # start working and block parent process until all processes done
        pool.close()
        pool.join()
        
        # complete the rest task
        self._collectResult(poolResults)
        self._writeResult()
        self._showStats()
    

    def _getResult(self, data):
        """(dict of JSON) -> (dict of lst, lst, lst)"""
        result = {}
        normalSamples = []
        tumorSamples = []
        
        for block in data:
            path = os.path.join(block["file_id"], block["file_name"])
            print("Current: " + path + "......")
            
            # e.g. 'TCGA-SX-A7SL-01A-11R-A355-07', '01A' is tumor, '11A' is normal
            submitterId = block['associated_entities'][0]['entity_submitter_id']
            if submitterId[13] == '0':
                tumorSamples.append(submitterId)
            else:
                normalSamples.append(submitterId)
                
            with gzip.open(path, 'rt') as f:
                lst = f.read().split('\n')
                if lst[-1] == '':
                    lst.pop()
                # e.g. {'ENSG00000185105.4' : [0.0059704985038]}
                for item in lst:
                    info = item.split('\t')
                    name = info[0]
                    value = info[1]
                    try:
                        result[name].append(value)
                    except KeyError:
                        result[name] = [value]
        return result, normalSamples, tumorSamples
    
    def _collectResult(self, poolResults):
        for result in poolResults:
        # get() return a list of tuple contain result
            valueDict, normals, tumors = result.get()[0]
                
            # combine result calculated from different process
            for key in valueDict:
                if key in self._results:
                    self._results[key] += valueDict[key]
                else:
                    self._results[key] = valueDict[key]
            
            self._normalSamples += normals
            self._tumorSamples += tumors
            
    def _writeResult(self):
        """Writing result to self._outFileName"""
        outFile = open(self._outFileName, "w")
        
        # Writing TGCA Normal first, then tumor
        if len(self._normalSamples) == 0:
            outFile.write("id")
        else:
            outFile.write('id\t' + '\t'.join(self._normalSamples))
        outFile.write('\t' + '\t'.join(self._tumorSamples) + '\n')
    
        # e.g. ENSG00000185105.4
        for key in self._results.keys():
            values = self._results[key]
            outFile.write(key + '\t')
            outFile.write('\t'.join(values) + '\n')
        outFile.close()
         
    def _showStats(self):
        print("\n---------------- Stats ----------------")
        print("Normal Count: " + str(len(self._normalSamples)))
        print("Tumor Count: " + str(len(self._tumorSamples)))
        print("------------------------------------------\n")
        print("Result has been generated to " + self.outFileName)
        print("\nProcess done, it takes: " +
              str(round(time.time() - self._startTime, 2)) + "sec\n\n"
              )
        print("Please wait, program exiting.......\n")
        
    
if __name__ == '__main__':
    analyzer = DataAnalyzer(sys.argv[1], "mRNA_Matrix.txt")
    analyzer.run()

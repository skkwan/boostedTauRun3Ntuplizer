To run the two file solution with AODSIM to find the parent (raw) and child (mini):
1. First place your analyzer in this folder and call it test-Analyzer.py
2. Next put the AODSIM files in inputFileListAODSIM.txt
3. Create the submit footer files using the commmand ```python dasqueryAODtoMiniRaw.py inputFileListAODSIM.txt```
For the next step PLEASE run a test job to see that it runs and outputs a file on /hdfs/ before submitting all jobs.
4. Create and submit the analyzers bash createAnalyzers.sh
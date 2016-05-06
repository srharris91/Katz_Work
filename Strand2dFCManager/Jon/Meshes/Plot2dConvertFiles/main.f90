PROGRAM main
CHARACTER(40)::inputFile
inputFile = "GenerateMesh.input"
WRITE(*,*)"OPENING ", Trim(AdjustL(inputFile))," as Input File"
call ReadPlot2Strand(inputFile)

WRITE(*,*) "PROGRAM COMPLETED!"

END PROGRAM main

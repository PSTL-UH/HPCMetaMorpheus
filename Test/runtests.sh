#!/bin/bash
if [ ! -f "elements.dat" ]; then
    ln -s ../mzlib/UsefulProteomicsDatabases/Datafiles/elements.dat .
fi

echo " "
echo "  BinGenerationTest"
make -s BinGenerationTest
./BinGenerationTest


echo " "
echo "  XLTest"
make -s XLTest
./XLTest

echo " "
echo "  XLSearchOutputTest"
make -s XLSearchOutputTest
./XLSearchOutputTest

echo " "
echo "  TestToml"
make -s TestToml
./TestToml

echo " "
echo "  IndexEngineTest"
make -s IndexEngineTest
./IndexEngineTest

echo " "
echo "  FdrTest"
make -s FdrTest
./FdrTest

#echo " "
#echo "  SearchEngineTests"
#make -s SearchEngineTests
#./SearchEngineTests

echo " "
echo "  CalibrationTests"
make -s CalibrationTests
./CalibrationTests


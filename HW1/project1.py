#
# Filename: project1.py
# Description: 
#
# Author: Gordon Kiesling
# VUNetID: Kiesligs
# Email: Gordon.s.kiesling@vanderbilt.edu
# Date: 1/29/16
#
# Honor Statement: "I pledge on my honor that I have neither given nor received unauthorized aid on this examination.”
#
# Caveats: 
#

from dna_strand import *


# Program main methoD
def main():
    # test some basic functions of the DNA_Strand class
    dna_1 = DNAStrand()  # start with two empty strands
    dna_2 = DNAStrand()

    # test length. This will automatically invoke the __len__() member function
    if len(dna_1) != 0:
        print("Length test failed for default constructor.")
    
    # or you can call  __len__() directly
    if dna_1.__len__() != 0:
        print("Length test failed for default constructor.")
    
    # test equality. this invokes __eq__()
    if not dna_1 == dna_2:
        print("Equality test failed.")
    
    # test should be reflective
    if not dna_2 == dna_1:
        print("Equality test failed.")
    
    # create a non-empty DNA
    dna = DNAStrand("ACTTGATTGGGTTGCTTGCC")

    print(dna.cleave("TTG"))

    # test length
    if len(dna) != 20:
        print("Length test failed for parameterized constructor")
    
    # test to string
    if dna.__str__() != "ACTTGATTGGGTTGCTTGCC":
        print("To string test failed.")
    
    # you can also use str() to automatically call __str__()
    # the following call is identical to the one above
    if "ACTTGATTGGGTTGCTTGCC" != str(dna):
        print("To string test failed.")
        
    # test indexing - invokes the __getitem__() method
    if dna[0] != "A":
        print("__getitem__() test failed.")
    try:
        dna[100]
        print("__getitem__() test failed.")
    except:
        pass
    
    # test indexing - __setitem__
    dna[0] = "X"
    if dna[0] != "X":
        print("__setitem__() test failed.")
    try:
        dna[100] = "X"
        print("__setitem__() test failed")
    except:
        pass

    # test search
    if dna.search("G") != 4:
        print("Search test failed.")
    if dna.search("G", 100) != -1:
        print("Search test failed.")

    # add your testing code here
    # be sure to test as thoroughly as possible
    
    # test in_range
    if not dna.in_range(1):
        print("__in_range__() test failed")
    if dna.in_range(100) != False:
        print("__in_range__() test failed")

    # test count_enzyme
    if dna.count_enzyme("TTG") != 4:
        print("count_enzyme test failed")

    # test append
    if dna.append(" is a dna strand, prof") != "ACTTGATTGGGTTGCTTGCC is a dna strand, prof":
        print("append test failed")

    # test cleave
    if dna.cleave("TTG") != -1:
        print("cleave test failed")

    # test cleave_all
    if dna != "ACTTGGGTTGCC":
        print("cleave_all test failed")

    # test splice
    if dna.splice("TTG,","bruh") != "ACTTGbruhGGTTGbruhCC":
        print("splice test failed")


    print("Done testing")


# The following allows the main method
# to run automatically when this file
# is passed to python
if __name__ == "__main__":
    main()

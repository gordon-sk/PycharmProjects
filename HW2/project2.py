"""
 File name: project1.py
 Author: Gordon Kiesling
 VUnetid: Kiesligs
 Email: Gordon.s.kiesling@vanderbilt.edu
 Class: CS2204
 Date: Monday, Feb 8 2016
 Honor statement: I have neither given nor received unauthorized aid on this work.
 Assignment Number: HW2
 Description: test driver for DNA_Strand.py
"""

from dna_strand import *

"""
Note: The way in which this assignment writes comments differs from the 
    first assignment. Python can make use of either system and you should
    be aware of both.
"""


def main():
    # test some basic functions of the DNA_Strand class
    dna_1 = DNAStrand()
    dna_2 = DNAStrand()
    
    if dna_1.size() != 0:
        print("Test 1 FAIL (size)")
    
    if not dna_1.__eq__(dna_2):
        print("Test 2 FAIL (eq)")
    
    if not (dna_1 == dna_2):   # == operator calls __eq__
        print("Test 2a FAIL (eq)")
    
    # test should be reflective
    if not dna_2.__eq__(dna_1):
        print("Test 3 FAIL (eq)")

    # create a non-empty DNA
    dna = DNAStrand("ABCCTG")
    
    # toString should return the contents as a string
    if str(dna) != "ABCCTG":
        print("Test 4 FAIL. (toString")

    # be sure to test as thoroughly as possible

    # ne should return true if different lengths or true if indexes are different
    dna_new1 = DNAStrand("ABCDEF")
    dna_new2 = DNAStrand("GHIJKL")
    if dna_new1.__ne__(dna_new2):
        print("Test 5 FAIL. (ne)")

    # testing set item
    dna.__setitem__(2,"Z")
    if dna[2] != "Z":
        print("Test 6 FAIL. (setitem)")

    # testing get item
    if dna.__getitem__(2) != "Z":
        print("Test 6 FAIL. (getitem)")

    # testing in_range
    if not dna.in_range(4):
        print("Test 7a FAIL (in_range true")
    if dna.in_range(100):
        print("Test 7b FAIL (in_range false)")

    # testing append
    append1 = DNAStrand("ABCDEF")
    append2 = DNAStrand("ABCDEF")
    append1.append(append2)
    if append1[6] != "A":
        print("Test 8 FAIL (append")

    # testing search
    if append1.search("A") != 0:
        print("Test 9a FAIL: (search whole string)")
    if append1.search("A",6) != 6:
        print("Test 9b FAIL: (search half string)")
    if append1.search("Z") != -1:
        print("Test 9c FAIL: (search for blank)")

    # testing count_enzyme
    countDNA = DNAStrand("ABCCTTGAHFOJVETTGFWJEFGCTTGATTGTTG")
    if countDNA.count("TTG") != 5:
        print("Test 10 FAIL: Count_enzyme")

    # testing cleave
    cleaveDNA = DNAStrand("AGCCTTGAGCATTGABTTG")
    if cleaveDNA.cleave("TTG") != ['A', 'G', 'C', 'C',
                                   'T', 'T', 'G', 'A',
                                   'B', 'T', 'T', 'G']:
        print("Test 11 FAIL: Cleave")

    # testing cleave_all
    cleaveallDNA = DNAStrand("ACTTGATTGGGTTGCTTGCC")
    if cleaveallDNA.cleave_all("TTG") != ['A', 'C', 'T', 'T', 'G', 'G',
                                          'G', 'T', 'T', 'G', 'C', 'C']:
        print("Test 11 FAIL : Cleave_all")

    # testing splice
    spliceDNA = DNAStrand("ACTTGATTGGGTTGCTTGCC")
    if spliceDNA.splice("TTG","bro") != ['A','C','T','T','G','b','r','o','G','G',
                                   'T','T','G','C','T','T','G','C','C']:
        print("Test 12 FAIL : Splice")



    print("Done testing")

if __name__ == '__main__':
    main()

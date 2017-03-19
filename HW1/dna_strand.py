#
# Filename: DNAStrand.py
# Description: This class defines DNA strings as objects and includes a variety
#              of methods that can both tell the user information about their
#              set of DNA and modify that set in similar ways that enzymes
#              found in nature do.
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


# Represents a DNA strand using strings to store the data
#
# To reduce errors in the starter files, return statements have been
# included beneath each method. Make the decision on whether the return
# statement is necessary for a particular method.
class DNAStrand(object):
    # Graduate students set to true
    GRAD_STUDENT = False

    # Constructor
    #
    #  @param (optional) 
    #    dna is a string DNA sequence
    #    defaults to empty string
    def __init__(self, dna=""):
        self._mDNA = dna

    # Returns the DNAStrand as a string
    #
    # @return
    def __str__(self):
        return str(self._mDNA)

    # Returns the length of the DNAStrand
    #
    # @return
    def __len__(self):
        return len(self._mDNA)

    # Compare this DNAStrand with rhs for equality.
    # @param rhs DNA to compare
    # @return
    #    True if the DNA strings are equal
    #    else, False
    def __eq__(self, rhs):
        return self._mDNA == rhs

    # Get an item in the DNAStrand at location index. Throws
    # IndexError if index is out of range, i.e., larger than the
    # current size() of the DNAStrand.
    def __getitem__(self, position):
        if -1 < position < len(self._mDNA):
            return self._mDNA[position]
        else:
            raise IndexError("Index out of range")

    # Set an item in the DNAStrand at location index. Throws
    # IndexError if index is out of range, i.e., larger than the
    # current size() of the DNAStrand.
    def __setitem__(self, position, letter):
        if -1 < position < len(self._mDNA):
            self._mDNA.replace(str(self._mDNA)[position], letter, 1)
            return
        else:
            raise IndexError("Index out of range")

    # Verifies that position is between 0 (inclusive) and the
    # length of the DNA sequence (exclusive)
    #
    # @param position
    # @return
    #    True if position within range
    #    False if not
    def in_range(self, position):
        return -1 < position < len(self._mDNA)

    # Locates the first instance of target
    #
    # @target character to locate
    # @position (optional) 
    #    the starting position
    #    defaults to 0
    # @return
    #    index of target if found
    #    -1 if target not found
    def search(self, target, position=0):
        return self._mDNA.find(target, position, len(self._mDNA))

    # Removes from current DNA strand the sequence between the *end* of the
    # first occurrence of passed target sequence (e.g. "TTG"), through the 
    # *end* of the second occurence of the target sequence.
    #
    # @pre e.g. "ACTTGACCTTGA" and target e.g. "TTG"
    # @post ACTTGA  (ACCTTG removed)
    # @param target target sequence
    # @param position (optional)
    #    the starting position
    #    defaults to 0
    # @return index following the first instance of target post-cleave 
    #    -1 if fewer than two instances of target after position
    #    (hint: helpful for cleaveAll)
    def cleave(self, target, position=0):
        length = len(target)
        start = self._mDNA.find(target)
        DNA = str(self._mDNA)

        part_1 = DNA[position:DNA.find(target)+target.__len__()]
        part_2 = DNA[DNA.find(target)+target.__len__():DNA.__len__()]
        part_3 = part_2[part_2.find(target)+target.__len__():DNA.__len__()]
        new_dna_string = part_1 + part_3


        return new_dna_string.find(target)

    # Removes from current DNA strand the sequence between pairs of target
    # sequence, i.e. from the end 1 through the end of 2, from the end of 3 
    # through the end of 4, etc, but NOT from the end of 2 through the end 3,
    # or from the end of 4 through the end of 5.
    # (Make sure that you understand the specification)
    #
    # 1-->2 and 3--->4
    #
    # @pre DNA e.g. ACTTGATTGGGTTGCTTGCC and target e.g. "TTG"
    # @post DNA ACTTGGGTTGCC (ATTG and CTTG removed)
    # @param target target sequence
    def cleave_all(self, target):
        DNA = str(self._mDNA)
        length = DNA.__len__()
        newDNA = ""
        split = (0,0,0)

        for x in range (0,length):
            if DNA.find(target)!= -1:
                split = DNA.partition(target)
                DNA = split[2][DNA.find(target)+(target.__len__() - 1):]
                newDNA = newDNA + split[0] + split[1]
        newDNA += DNA

        self._mDNA = newDNA
        return

    # Counts the number of occurrences of a single character target sequence
    # in the DNA strand
    #
    # @param target single enzyme as character
    # @return number of occurrences of target
    def count_enzyme(self, target):
        DNA = str(self._mDNA)
        return int(DNA.count(target))

    # Append the characters of the parameter to the end of the current DNA
    # @pre DNA e.g. ACTTGA and rhs e.g. ACCTG
    # @post DNA ACTTGAACCTG
    # @param rhs the DNA strand to append to the right hand side of existing
    def append(self, rhs):
        DNA = str(self._mDNA)
        return str(DNA + rhs)

    # Replaces the DNA sequence between the end of the first occurrence of
    # target and the end of the second occurrence of target with a new
    # DNA sequence.
    # @post
    #    if two targets found, target is replaced by replacement
    #    else, no change
    # @param target specifies the DNA sequence to replace
    # @param replacement DNA sequence replacing the removed sequence
    # @param position (optional)
    #    position to begin searching for target
    #    defaults to 0
    # @return 
    #    index following the first instance of replacement post-splice
    #    -1 if fewer than two instances of target after position
    def splice(self, target, replacement, position=0):

        length = len(target)
        start = self._mDNA.find(target)
        DNA = str(self._mDNA)

        part_1 = DNA[position:DNA.find(target)+target.__len__()]
        part_2 = DNA[DNA.find(target)+target.__len__():DNA.__len__()]
        part_3 = part_2[part_2.find(target)+target.__len__():DNA.__len__()]
        new_dna_string = part_1 + part_3 + replacement
        self._mDNA = new_dna_string

        if self._mDNA.count(target) < 2:
            return -1
        else:
            return self._mDNA[self._mDNA.find(replacement):]

    #########################
    # THE FOLLOWING METHODS ARE FOR GRAD STUDENTS ONLY
    if GRAD_STUDENT:
        MARKER = "#"  # you can reference this variable by using: self.__class__.MARKER

        # Find all non-overlapping occurrences of the target sequence and insert
        # the '#' marker AFTER each of them.
        #
        # @pre DNA e.g. ACTTGATTGGGTTGCTTGCC and target e.g. "TTG"
        # @post DNA ACTTG#ATTG#GGTTG#CTTG#CC
        def insert_marker(self, target):
            # YOUR CODE GOES HERE
            return

        # Delete all the '#' markers that occur in the strand
        # @pre DNA e.g. ACTTG#ATTG#GGTTG#CTTG#CC
        # @post DNA ACTTGATTGGGTTGCTTGCC
        def delete_marker(self):
            # YOUR CODE GOES HERE
            return

        # Returns the number of markers in the strand
        # @return
        def count_marker(self):
            # YOUR CODE GOES HERE
            return

        # Similar to cleaveAll finds pairs of targets in current DNA strand and replaces
        # the sequence between the end of the first target through the end of the 
        # second with the insertSequence. If two instances of the target are not found, 
        # then no changes are made.
        def splice_all(self, target, replacement):
            # YOUR CODE GOES HERE
            return

            # endif GRAD_STUDENT
            # END GRAD STUDENT ONLY METHODS
            #########################

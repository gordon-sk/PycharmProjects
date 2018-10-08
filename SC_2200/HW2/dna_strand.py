"""
 File name: DNAStrand.py
 Author: Gordon Kiesling
 VUnetid: Kiesligs
 Email: Gordon.s.kiesling@vanderbilt.edu
 Class: CS2204
 Date: Monday, Feb 8 2016
 Honor statement: I have neither given nor recieved unauthorized aid on this work.
 Assignment Number: HW2
 Description: This will be a DNAStrand implemented with a python list.
"""


class DNAStrand(object):
    """
    This class represents a DNA strand with an internal list of characters.
    Supported operations include searching, cleaving and splicing.
    """
    # NOTE: change this constant to True if you are a graduate student
    GRAD_STUDENT = False

    def __init__(self, data=None):
        """
        Constructor
        param data accepts an existing sequence of characters (string, list, tuple, or array)
        """
        # initialize the data list
        # this method is done for you
        if data is None:
            self._myDNA = []
        else:
            self._myDNA = [ch for ch in data]
            
    def __str__(self):
        """
        return string equivalent of the DNA
        """
        # your code here
        strdna = ""
        for x in range(self._myDNA.__len__()):
            strdna = strdna + self._myDNA[x]
        return strdna
    
    def __len__(self):
        """
        return the size of the DNAStrand
        """
        # your code here
        return self._myDNA.__len__()
    
    # alternate definition for size accessor
    # don't modify this line
    size = __len__
    
    def __iter__(self):
        """
        return an iterator for this DNAStrand
        """
        return self._myDNA.__iter__()    # this one is done for you
    
    def __eq__(self, rhs):
        """
        __eq__
        Compare this DNAStrand with rhs for equality. Returns true if the
        size()'s of the two DNAStrands are equal and all the elements of the
        list are equal, else false.
        """
        eq1 = self._myDNA.__len__() == rhs.__len__()
        eq2 = True
        if eq1:
            for x in range(self._myDNA.__len__()):
                if self._myDNA[x] != rhs[x]:
                    eq2 = False
                else:
                    eq2 = True
        return eq1 & eq2

    def __ne__(self, rhs):
        """
        Not equals operator.
        returns the negation of the equality operator for
        two given objects
        """
        value = 0
        if  self._myDNA.__eq__(rhs):
            value = False
        else:
            value = True

        return value
        
    def __setitem__(self, index, letter):
        """
        Set a letter in the DNAStrand at location index. Throws
        IndexError if index is out of range, i.e., larger than the
        current size() of the DNAStrand or less than zero. Throws
        ValueError if the letter is not a single character.
        """
        if index > self._myDNA.__len__():
            raise IndexError("Index out of range, chief!")
        if str(index).__len__() != 1:
            raise ValueError("Index is larger than single character")

        self._myDNA[index] = letter
        return

    # alternate definition for item mutator
    # don't modify this line
    set = __setitem__
    
    def __getitem__(self, index):
        """
        Get a letter in the DNAStrand at location index. Throws
        IndexError if index is out of range, i.e., larger than the
        current size() of the DNAStrand or less than zero.
        """
        if not 0 < index < self._myDNA.__len__():
            raise IndexError("Index out of range, chief!")
        return self._myDNA[index]
    
    # alternate definition for item accessor
    # don't modify this line
    get = __getitem__
    
    def in_range(self, position):
        """
        inRange : helper function
        Returns true if position is within range, i.e., 0 <= position < length of DNA strand
        else returns false.
        :param position: location in strand
        """
        return 0 < position < self._myDNA.__len__()
        
    def search(self, target, position=0):
        """
        search
        search with start position specified
        Look for target in current DNA strand and return index.
        Return -1 if not found.
        Pre: target is a string
        :param position: default 0, can be made larger to start search
                         farther intro strand
        :param target: the character being searched for
        """
        try:
            value = self._myDNA.index(target, position)
        except ValueError:
            value = -1
        return value

    # Alternate definition of search
    # more consistent with python collection semantics
    # don't modify this line
    index = search
    
    def cleave(self, target, position=0):
        """
        cleave
        Removes from current DNA strand the sequence between the end of the 
        first occurrence of passed target sequence (e.g. "TTG"), through the end
        of the second occurrence of the target sequence.
        Returns the first index after the cleave
        pre: Array e.g. ACTTGACCTTGA and target is a string e.g. "TTG"
        post: ACTTGA  (ACCTTG removed), index 5 returned
        :param position: default 0, can be upped to start search later in string
        :param target: DNA mini-strand we are looking for
        """
        if target.__len__() == 0:
            raise ValueError("Target is an empty string!")
        record = []
        x = position
        while x < self._myDNA.__len__():
            if self._myDNA[x] == target[0]:
                record.append(x)
                x = x + target.__len__()
            else:
                x += 1
        return self._myDNA[position:record[0]+target.__len__()] + self._myDNA[record[1]+target.__len__():]

    def cleave_all(self, target):
        """
        cleaveAll
        Removes from current DNA strand the sequence between pairs of target 
        sequence, i.e. from the end 1 through the end of 2, from the end of 3 
        through the end of 4, etc, but NOT from the end of 2 through the end 3,
        or from the end of 4 through the end of 5.
        (Make sure that you understand the specification)
        pre: Array e.g. ACTTGATTGGGTTGCTTGCC and target (a string) e.g. "TTG"
        post: ACTTGGGTTGCC (ATTG and CTTG removed)
        :param target: DNA string we are searching for
        """
        if target.__len__() == 0:
            raise ValueError("Target is an empty string!")
        record = []
        x = 0
        while x < self._myDNA.__len__():
            if self._myDNA[x] == target[0]:
                record.append(x)
                x = x + target.__len__()
            else:
                x += 1
        newdna = self._myDNA[:record[0]+target.__len__()]
        y = 1
        while y < record.__len__()-1:
            newdna = newdna + self._myDNA[record[y]+target.__len__():record[y+1]+target.__len__()]
            y += 2
        newdna = newdna + self._myDNA[record[y]+target.__len__():]

        return newdna
        
    def count_enzyme(self, target):
        """
        countEnzyme
        Counts the number of non-overlapping occurrences of a target sequence
        in the DNA strand
        Eg, the target "AAA" appears 3 non-overlapping times in the DNA "AAAAAAAAAA"
        Pre: target is a nonempty string. Raise ValueError if target length is zero
        :param target: DNA target we are searching for
        """
        if target.__len__() == 0:
            raise ValueError("Target is an empty string!")
        record = []
        x = 0
        while x < self._myDNA.__len__():
            if self._myDNA[x] == target[0]:
                record.append(x)
                x = x + target.__len__()
            else:
                x += 1
        count = 0
        match = False
        for y in range(record.__len__()):
            match = False
            for z in range(target.__len__()):
                if self._myDNA[record[y] + z] == target[z]:
                    match = True
                else:
                    match = False
                    break
            if match:
                count += 1
        return count

    # Alternate definition of count
    # more consistant with python collection semantics
    # don't modify this line
    count = count_enzyme
    
    def append(self, rhs):
        """
        Append the characters of the parameter to the end of the current DNA.
        Example: if _myDNA contained ACTTGA and 'ACCTG' was received as a parameter, 
        then afterward _myDNA will contain ACTTGAACCTG
        Pre: rhs is a string parameter
        :param rhs: list we are appending
        """
        self._myDNA = self._myDNA + [ch for ch in rhs]
        return

    def splice(self, target, replacement, position=0):
        """
        splice (accepts 2 Strings representing sequences)
        finds first pair of targets in current DNA strand and replaces
        the sequence between the end of the first target through the end of the 
        second with the replacement. 
        If two instances of the target are not found, 
        then no changes are made.
        Returns the index of the position after the splice, or -1 if no splice was performed
        Pre: target & replacement are strings
        :param position: default 0, can be upped to start search farther intro list
        :param replacement: string that will be implimented into list
        :param target: DNA string we are searching for
        """
        if target.__len__() == 0:
            raise ValueError("Target is an empty string!")
        record = []
        x = 0
        while x < self._myDNA.__len__():
            if target[0] == self._myDNA[x]:
                record.append(x)
                x = x + target.__len__()
            else:
                x += 1
        return self._myDNA[position:record[0]+target.__len__()] + [ch for ch in replacement] + self._myDNA[record[1]+
                                                                                                   target.__len__():]

    # ########################/
    # THE FOLLOWING METHODS ARE FOR GRAD STUDENTS ONLY
    if GRAD_STUDENT:
        _MARKER = '#'
        
        def insert_marker(self, target):
            """
            insertMarker
            Find all non-overlapping occurrences of the target sequence and insert
            the '#' marker AFTER each of them.
            pre: Array e.g. ACTTGATTGGGTTGCTTGCC and target e.g. "TTG"
            post: ACTTG#ATTG#GGTTG#CTTG#CC
            """
            # your code here
            pass
            
        def delete_marker(self):
            """
            deleteMarker
            Delete all the '#' markers that occur in the strand, shifting data to the
            left as appropriate
            """
            # your code here
            pass
            
        def count_marker(self):
            """
            count_MARKER
            return the number of _MARKERs that exist in the strand
            """
            # your code here
            return -1
            
        def splice_all(self, target, replacement):
            """
            spliceAll (accepts 2 Strings representing sequences)
            Similar to cleaveAll finds pairs of targets in current DNA strand and replaces
            the sequence between the end of the first target through the end of the 
            second with the insertSequence. If two instances of the target are not found, 
            then no changes are made.
            """
            # your code here
                
    # END GRAD STUDENT ONLY METHODS
    #########################/
# File name: Gradebook.py
# Author: Gordon Kiesling
# VUnetid: Kiesligs
# Email: Gordon.s.kiesling@gvanderbilt.edu
# Class: CS2204
# Assignment Number: 7
# Honor statement: I have neither given nor received unauthorized aid on this work
# Description: Read in a gradesheet in csv format and compute each student's
#   weighted average, class rank, letter grade.

import csv
import math

# Static weight definitions
proj_pct = 0.45  # weight given homework assignments
exam1_pct = 0.10  # weight given exam #1
exam2_pct = 0.10  # weight given exam #2
exam3_pct = 0.15  # weight given exam #3
exam4_pct = 0.20  # weight given the final exam


class Gradebook:
    # Reads data from the provided .csv file and initializes an 
    # internal 2d list to store the values.
    # Perform required pre-processing of raw data:
    #  1) any row that is longer than the first row is truncated (a warning message is also printed)
    #  2) any row that is shorter than the first row is padded with '0'
    #  3) all empty strings or single blank strings in the grades are replaced with '0'
    def __init__(self, filename):
        self.gradesheet = []
        file = open(filename)
        data = file.readline()
        x = 0
        while data is not "":
            self.gradesheet.append(data.split(','))
            data = file.readline()
            x += 1
        self.clean_up(self.gradesheet)  # Our three post-conditions in separate helper methods
        file.close()
        pass

    # Ensures all lines in the gradesheet list are no longer than the first line
    # Replaces all empty strings with the value 0
    # Ensures shorter lines are extended with 0s
    # pre: gradesheet is a list of lists representing rows in CSV file of grades
    # post: all rows are made to length of first row, unacceptable values converted to 0
    def clean_up(self, gradesheet):
        max_len = gradesheet[0].__len__()  # setting the target length

        for x in range(gradesheet.__len__()):  # shortening / lengthening
            if gradesheet[x].__len__() > max_len:
                new_row = []
                for y in range(max_len):
                    new_row.append(gradesheet[x][y])
                gradesheet[x] = new_row
                print("Error: Line " + repr(x) + " in grades longer than expected, has been truncated\n")
            elif gradesheet[x].__len__() < max_len:
                for y in range(gradesheet[x].__len__(), max_len):
                    gradesheet[x].append(0)

            for y in range(gradesheet[x].__len__()):  # checking for empty values
                if gradesheet[x][y] == '' or gradesheet[x][y] == ' ':
                    gradesheet[x][y] = 0

            for z in range(self.gradesheet.__len__()):
                gradesheet[z][28] = str(gradesheet[z][28]).rstrip()

    # pre: internal gradebook has been created/initialized
    # post: weighted average, grade, and class rank have been computed for all students
    # The order of method calls MUST be maintained
    def compute_grades(self):
        self.compute_weighted_average()
        self.compute_class_rank()
        self.assign_grades()

    # pre:  internal gradebook has been created/initialized
    # post: Internal gradebook now contains an
    #       additional column containing each 
    #       student's weighted average. This column
    #       has a 'Weighted Average' header.
    # note: Some students fail to submit assignments,
    #       ensure the these students receive 0s for
    #       those assignments.
    # return: list of weighted averages
    def compute_weighted_average(self):
        self.avg_list = []
        proj_pt_pos = 0
        for x in range(3, 25):
            proj_pt_pos += int(self.gradesheet[1][x])
        for x in range(2, self.gradesheet.__len__()):
            proj_pt_earned = 0
            for y in range(3, 25):
                proj_pt_earned += int(self.gradesheet[x][y])
            self.avg_list.append(100 * ((proj_pt_earned / proj_pt_pos) * proj_pct +
                                        (int(self.gradesheet[x][25]) / 100) * exam1_pct +
                                        (int(self.gradesheet[x][26]) / 100) * exam2_pct +
                                        (int(self.gradesheet[x][27]) / 100) * exam3_pct +
                                        (int(self.gradesheet[x][28]) / 100) * exam4_pct))

        for x in range(self.avg_list.__len__()):  # rounding out numbers
            self.avg_list[x] = round(self.avg_list[x], 2)
        for x in range(self.avg_list.__len__()):  # appending to gradesheet
            self.gradesheet[x + 2].append(self.avg_list[x])
        return self.avg_list

    # pre:  Requires that the weighted averages have
    #       been calculated.
    # post: Internal gradebook now contains an
    #       additional column containing each
    #       student's class rank.
    # return: list of ranks
    def compute_class_rank(self):
        temp_avg = self.avg_list
        temp_avg.sort(reverse=True)
        self.rank_list = [0 for num in range(temp_avg.__len__())]
        for x in range(temp_avg.__len__()):
            for y in range(temp_avg.__len__()):
                if self.gradesheet[x + 2][29] == temp_avg[y]:
                    self.rank_list[x] = y + 1
                    self.gradesheet[x + 2].append(y + 1)
        return self.rank_list

    # pre:  Requires that the weighted averages and
    #       class ranks have been calculated
    # post: Internal gradebook now contains an
    #       additional column containing each
    #       student's letter grade.
    # return: list of grades
    def assign_grades(self):
        self.grades_list = [0 for num in range(self.gradesheet.__len__() - 2)]
        for x in range(2, self.gradesheet.__len__()):
            if math.ceil(self.gradesheet[x][29]) >= 90:  # tedious code ahead...
                self.gradesheet[x].append("A")  # there was no other way
                self.grades_list[x - 2] = "A"
            elif 80 <= math.ceil(self.gradesheet[x][29]) < 90:
                self.gradesheet[x].append("B")
                self.grades_list[x - 2] = "B"
            elif 70 <= math.ceil(self.gradesheet[x][29]) < 80:
                self.gradesheet[x].append("C")
                self.grades_list[x - 2] = "C"
            elif 60 <= math.ceil(self.gradesheet[x][29]) < 70:
                self.gradesheet[x].append("D")
                self.grades_list[x - 2] = "D"
            else:
                self.gradesheet[x].append("F")
                self.grades_list[x - 2] = "F"
        return self.grades_list

    # pre: none
    # post: Gradesheet has been saved to csv file
    # param: Output filename
    def write_grades(self, filename):
        file = open(filename, 'w')
        for x in range(self.gradesheet.__len__()):
            for y in range(self.gradesheet[x].__len__()):
                if y is not 0:
                    file.write(",")
                file.write(str(self.gradesheet[x][y]))
            file.write("\n")
        file.close()

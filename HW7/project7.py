# Driver script for project 7

from Gradebook import Gradebook


def main():
    # Create a Gradebook object from the specified file
    gb = Gradebook('grades.csv')
    
    # Compute the grades
    gb.compute_grades()
    
    # Write the grades (with new columns) to an updated file
    gb.write_grades('updated_grades.csv')

    # For quick progress checking
    print('Wt.Avg  Rank Grade')
    for row in gb.gradesheet[2:]:
        print('%05.2f%%%5d%5s' % (row[-3], row[-2], row[-1]))


if __name__ == '__main__':
    main()

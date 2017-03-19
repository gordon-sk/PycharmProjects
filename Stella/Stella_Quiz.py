"""
Name: Gordon Kiesling
File: Easy to take electronic form of stella survey by USNO
"""

def quiz():
    file_name = "Stella Users Survey Results May 2016 USS Donald Cook"
    file = open(file_name, "a")

    print()
    file.write("\n\n" + input("Please enter your rate and last name (Please do not include either slash character): ") + "\n\n")

    print("Hey, these are a few questions from the Naval Observatory regarding how often and for what")
    print("your ship is using Stella for. They're all brief and none are open ended -- should take around 5 minutes. \n")

    file.write("Ship: " + input("What is the name and type of your current ship?: ") + "\n")
    file.write("how often are cel fixes taken: " + input("How often are celestial fixes obtained on your ship, if weather conditions permit?: ") + "\n")
    file.write("Do you enter all of your ship's course and speed changes into Stella?: " + input("Do you enter all of your ship's course and speed changes into Stella? (Y/N): ") + "\n")
    file.write("Do you record/input weather data when doing a sight reduction?: " + input("Do you record/input weather data when doing a sight reduction? (Y/N): ") + "\n")
    file.write("What version of Windows is used to run Stella?: "+ input("What version of Windows is used to run Stella? XP, Vista, Win7 or other?: ")  + "\n")
    file.write("Does more than one person access Stella on the same Windows account? (Y/N): " + input("Does more than one person access Stella on the same Windows account? (Y/N): ")  + "\n")
    file.write("What is the highest privilege level of users? Admin, standard, or other?: " + input("What is the highest privilege level of users? Admin, standard, or other?: ")  + "\n")
    file.write("What works best for obtaining software and updates? CDs, web, or other?: " + input("What works best for obtaining software and updates? CDs, web, or other?: ")  + "\n")
    print("\nFunctions: Please indicate with the numbers 1-5 how useful each of the \nfollowing tasks are, with 1 = not used and 5 = essential for navigation\n")
    print("Almanac")
    file.write("Almanac\n1. Sun / Moon / Aries: " + input("1. Sun / Moon / Aries: ") + "\n")
    file.write("2. Planets: " + input("2. Planets: ")  + "\n")
    file.write("3. Stars: " + input("3. Stars: ")  + "\n")
    file.write("4. All: " + input("4. All: ")  + "\n")
    print("Position Update")
    file.write("Position Update\n1. Update DR: " + input("1. Update DR: ")  + "\n")
    file.write("2. External Fix: " + input("2. External Fix: ")  + "\n")
    file.write("3. Waypoint: " + input("3. Waypoint: ")  + "\n")
    print("Rise, Set & Transit\n")
    file.write("Rise, Set & Transit\n1. Fixed Site: " + input("1. Fixed Site: ")  + "\n")
    file.write("2. Underway: " + input("2. Underway: ")  + "\n")
    print("Sight Planning")
    file.write("Sight Planning\n" + input("1. Sky Chart: ")  + "\n")
    file.write("2. Selected Stars: " + input("2. Selected Stars: ")  + "\n")
    print("Gyro / Compass Error")
    file.write("Gyro / Compass Error\n1. Planning: " + input("1. Planning: ")  + "\n")
    file.write("2. Determination: " + input("2. Determination: ")  + "\n")
    print("Sight Reduction ")
    file.write("Sight Reduction\n1. Record Observations: " + input("1. Record Observations: ")  + "\n")
    file.write("2. Compute Fix: " + input("2. Compute Fix: ")  + "\n")
    print("Edit")
    file.write("Edit\n1. Correction Observation: " + input("1. Correct Observation: ")  + "\n")
    file.write("2. Correct Reference Position:" + input("2. Correct Reference Position:")  + "\n")
    file.write("3. Delete Observation: " + input("3. Delete Observation: ")  + "\n")
    file.write("4. Delete Reference Position: " + input("4. Delete Reference Position: ")  + "\n")
    file.write("5. Add Comment: " + input("5. Add Comment: ")  + "\n")
    print("Windows")
    file.write("Windows\n1.Worksheet: " + input("1. Worksheet: ")  + "\n")
    file.write("2. Current Log: " + input("2. Current Log: ")  + "\n")
    file.write("3. Current Track: " + input("3. Current Track: ")  + "\n")
    file.write("4. Sight Plan Plot: " + input("4. Sight Plan Plot: ")  + "\n")
    file.write("5. Sight Reduction Plot: " + input("5. Sight Reduction Plot: ")  + "\n")
    file.write("6. Strip Form: " + input("6. Strip Form: ")  + "\n")

    print("\n\nThat's all, thank you.")
    file.close()

if __name__ == "__main__":
    quiz()
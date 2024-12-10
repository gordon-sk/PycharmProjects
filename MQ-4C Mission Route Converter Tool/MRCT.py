# This script was created by LT Gordon "Yorkie" Kiesling in late 2024.
# This script was validated by XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX (give yourself some credit here guys!)
# The purpose of this script is to provide MQ-4C Air Vehicle Pilots the ability to quickly convert a mission plan Excel
# document (required for flight by NATOPS) into a .kml file. These files are broadly compatible with a wide range of
# applications, but the specific use case I had in mind when I created this script was for pilots to upload them into
# AFWEBS on SIPR so that they could more effectively utilize the weather overlays available there.
# It is intended to be run on Python 3.12, and requires only native libraries to that Python distro.
# If problems appear, I can probably be reached at gordon.s.kiesling.mil@us.navy.mil
# If you're reading this on an MCS, edit this file at your own risk.

import xml.etree.ElementTree as Element_Tree  # for constructing the .kml file
import os  # for making the terminal window pretty colors
import csv  # for parsing the .csv data the user will provide
import re  # For input sanitation
import sys  # for buffer clearing to streamline user input
from datetime import datetime, timezone  # for timestamping the output

# TODO: Pretty colors : )
# TODO: Maybe give the user the option to delete every csv and kml in the folder, if it's gotten crowded
# TODO: Add working area coordinates

RED = "\033[31m"
YELLOW = "\033[33m"
BLACK_BG = "\033[40m"  # Black background
BRIGHT_GREEN_TEXT = "\033[92m"  # Bright green text
RESET = "\033[0m"  # Resets color back to default
os.system('color 0A')
try:
    os.system('mode con: cols=175 lines=50')
except Exception as e:
    print(f"Error resizing: {e}")

# Gets the file name
# Opens the file
# Pulls the relevant coordinates
# Gets the user's current posit (optional) and inserts it into the coordinates list
# Creates a kml out of the coordinates
def main():

    csv_file_name = get_file()
    with open(csv_file_name, mode='r', newline='', encoding='utf-8') as file:
        data = csv.reader(file)
        coordinates, route_name = get_user_route(data)

    if not route_name == "NOT FOUND":
        current_pos = get_current_pos()
        if current_pos is not None:
            coordinates["UA POS"] = current_pos
        create_kml(coordinates, route_name)

    sys.stdin.flush()  # clears the buffer -- lingering \n problems caused user to have to sometimes hit enter twice.
    input("Press enter to close window.")


# FUNCTION get_file()
# Finds the directory the script is stored in, scans it for .csvs, and asks the user which one they'd like to work
# with. If no files are found, the program exits. Performs input sanitation on user's choice.
# Returns a string object representing the chosen file's name.
def get_file():
    script_dir = os.path.dirname(os.path.abspath(__file__))  # Finds the directory where this script is stored
    os.chdir(script_dir)  # Sets the working directory to script location
    print("This tool is stored in: " + os.getcwd())
    print("Please make sure you save the .csv file to that folder!\n")

    file_list = os.listdir(os.getcwd())  # Getting all the files in the directory
    csv_dict = {}
    file_choice_iterator = 1
    print("Here is a list of .csv files stored with this tool:")
    for file in file_list:
        file_name = file.split(".")
        # filtering for .csvs
        if file_name[len(file_name) - 1] == "csv":
            csv_dict[file_choice_iterator] = file
            file_choice_iterator += 1
    # if we don't find any...
    if len(csv_dict) == 0:
        print("Error: No .csv files found.")
        print("Please close this window and try again.")
    else:
        for i in range(len(csv_dict)):
            print(str(i + 1) + ": " + csv_dict[i + 1])

        while True:
            choice = input("Which file would you like to create a KML from? Type the number and hit enter: ")
            try:
                choice_int = int(choice)  # Convert input to an integer
                # We make sure the user chose a valid number
                if choice_int > len(csv_dict) or choice_int <= 0:
                    print(
                        f"Invalid input -- selection out of range. Please choose a number between 1 "
                        f"and {len(csv_dict)}.\n")
                else:
                    selected_file = csv_dict[choice_int]
                    print("File selected: " + selected_file + "\n")
                    return selected_file
            # Input has to be a number in the first place
            except ValueError:
                print("Invalid input -- must be a number. Please try again.\n")


# TODO: MTD copy paste support
# FUNCTION get_current_pos()
# Asks the user if they'd like to enter their current position, either manually or by copy and paste from the MTD.
# Use cases: User is far stitched waypoint, user is maneuvering around weather
# Skip cases: User is on deck, user position coincident with stitched waypoint
# Performs input sanitation and returns the current position as a .kml compliant tuple object
# with altitude set to 0.
#
# Example copy/paste: location;33.8936117994102;31.1578366687186
def get_current_pos():
    print("\nOptional: Enter your current lat/long. Press enter to skip.")
    print("You can right click on the MTD, hover over \"Tracks\", and click \"Copy location\", then paste it here with "
          "CTRL + V.")
    print("Alternatively, use the format on your PFD: XX:XX.XXN YYY:YY.YYE")
    pattern_1 = r"^\d{2}:\d{2}\.\d{2}[NS] \d{2,3}:\d{2}\.\d{2}[EW]$"
    pattern_2 = r"^LOCATION;-?[0-9]+\.[0-9]+;-?[0-9]+\.[0-9]+$"

    while True:
        sys.stdin.flush()  # clears the buffer -- lingering \n problems caused user to have to sometimes hit enter twice
        current_pos = input("\nCurrent Pos: ").strip().upper()
        # If the user skips this step:
        if current_pos == "":
            return None
        else:
            # manual entry option
            if re.match(pattern_1, current_pos):
                print("Position recorded: " + current_pos.upper())
                # we begin to turn it from a string into a tuple object compatible with the .kml creation process
                lat, lon = current_pos.split(" ")
                # we use regex to split the different parts of the latitude
                lat_parts = re.match(r"(\d+):(\d+)\.(\d+)([NS])", lat)
                # converting from strings to floats
                degrees = float(lat_parts.group(1))
                minutes = float(lat_parts.group(2))
                seconds = float(lat_parts.group(3))
                # checking if we're south of the equator
                direction = lat_parts.group(4)
                # combining into pure decimal form -- .kml compatible format is (-)XX.XXXX
                lat_decimal = degrees + minutes / 60 + seconds / 3600
                # Making it negative if relevant
                if direction in "S":
                    lat_decimal *= -1

                lon_parts = re.match(r"(\d+):(\d+)\.(\d+)([WE])", lon)
                degrees = float(lon_parts.group(1))
                minutes = float(lon_parts.group(2))
                seconds = float(lon_parts.group(3))
                direction = lon_parts.group(4)
                # combining into pure decimal form -- .kml compatible format is (-)YYY.YYYY
                lon_decimal = degrees + minutes / 60 + seconds / 3600
                if direction in "W":
                    lon_decimal *= -1

                # Declare and return the tuple
                current_pos_tuple = (lat_decimal, lon_decimal, 0)
                return current_pos_tuple

            # copy and paste option
            elif re.match(pattern_2, current_pos):
                lat = float(current_pos.split(";")[1])
                lon = float(current_pos.split(";")[2])
                current_pos_tuple = (lat, lon, 0)
                return current_pos_tuple
            else:
                print("Incorrect format detected. Please try again.\n")
                # Restarts, due to while True: logic above


# TODO: Optional C4 (go-around) data integration
# FUNCTION get_user_route()
# This function asks the user what route they would like to map, starting with a waypoint they specify.
# It does basic input sanitation, then parses through the .csv and pulls the relevant coordinates.
# INPUT data: a csv.reader object parseable by pure Python. Very handy, as you can read below.
# OUTPUT coord_list: a list object containing coordinates in the tuple format (lat, lon, altitude)
# OUTPUT route_name: a string object associated with the route the user specified.
def get_user_route(data):
    print("What route would you like to map? Enter the waypoint you intend to stitch to."
          "The .kml output will run from that waypoint to on-deck.")
    print("Include the leading character: P, A, L, R, Z, G. This is not case sensitive.")

    pattern = r'^[PALRZG]\d{1,5}$'
    route_name = ""
    coord_dict = {"UA POS": (0, 0, 0)}

    while True:
        starting_waypoint = input("Initial route waypoint: ").upper()
        print()
        if re.match(pattern, starting_waypoint):
            waypoint_found = False
            for row in data:

                if row[0] == "Route Name:":
                    route_name = row[3]

                if row[0] == starting_waypoint:
                    print("Starting waypoint found")
                    waypoint_found = True

                if waypoint_found:
                    waypoint_comment = row[2].split("\n")
                    coordinates = row[3].split('\n')
                    coordinates[0] = coordinates[0].strip().upper()
                    coordinates[2] = coordinates[2].strip().upper()

                    # Pulling the altitude data... Kind of useless for now but why not?
                    altitude = 0
                    if not len(row[7].split('\n')) == 1:
                        altitude = row[7].split('\n')[-1]

                    # Regex patterns
                    coord_format_1 = r"[NS] \d{1,3} \d{2} \d{2}\.\d{2}"  # ex:     N 34 06 58.80
                    coord_format_2 = r"[NS] \d{1,3}\.\d+"  # ex:     N 34.11633
                    coord_format_3 = r"[NS] \d{1,3} \d{2}\.\d+"  # ex:     N 26 20.5466

                    # KML coordinate format must comply with the following example:
                    # coordinates = [
                    #     (34.11633, -118.35258),  # (latitude, longitude)
                    #     (34.04556, -118.25078),
                    #     (34.10123, -118.34567)
                    # ]

                    split_lat = coordinates[0].split(" ")
                    split_long = coordinates[2].split(" ")

                    if re.match(coord_format_1, coordinates[0]):

                        lat = float(split_lat[1]) + float(split_lat[2]) / 60 + float(split_lat[3]) / 3600
                        long = float(split_long[1]) + float(split_long[2]) / 60 + float(split_long[3]) / 3600
                        if split_lat[0] == 'S':
                            lat *= -1
                        if split_long[0] == 'W':
                            long *= -1
                        coord_dict[row[0]] = (lat, long, altitude)

                    elif re.match(coord_format_2, coordinates[0]):

                        lat = float(split_lat[1])
                        long = float(split_long[1])
                        if split_lat[0] == 'S':
                            lat *= -1
                        if split_long[0] == 'W':
                            long *= -1
                        coord_dict[row[0]] = (lat, long, altitude)

                    elif re.match(coord_format_3, coordinates[0]):

                        lat = float(split_lat[1]) + float(split_lat[2]) / 60
                        long = float(split_long[1]) + float(split_long[2]) / 60
                        if split_lat[0] == 'S':
                            lat *= -1
                        if split_long[0] == 'W':
                            long *= -1
                        coord_dict[row[0]] = (lat, long, altitude)

                    elif coordinates[0].upper() == "LATITUDE":
                        print("Skipping blank row...")
                    else:
                        raise ValueError("Unknown coordinate formatting detected -- please let Yorkie know what mission"
                                         " plan you are using at gordon.s.kiesling.mil@us.navy.mil")

                    if waypoint_comment[len(waypoint_comment) - 1] == "Touchdown":
                        print("Touchdown point found: " + row[0] + ".")
                        break

            if len(coord_dict) > 1:
                return coord_dict, route_name
            else:
                print("Waypoint not found. Please close this window and try again.")
                return coord_dict, "NOT FOUND"

        else:
            print("Incorrect format detected. ")
            # Restarts, due to while True: logic above


# FUNCTION create_kml(coordinates, output_label)
# Function to create KML with lines connecting the coordinates provided.
# Icons are rendered invisible, as I feel they're more distraction than help.
# Function takes the following inputs:
# INPUT coordinates: a list object containing tuples representing coordinates the aircraft will fly-by.
# INPUT output_label: a string representing the mission route specified by the user. The mission route is the "header"
# for the specific subsection of the overall mission plan the aircraft has loaded.
# Outputs: None. Once the file is created, the user can close the window.
def create_kml(coordinates, output_label):
    # Create the root element <kml> and its namespace
    kml = Element_Tree.Element("kml", xmlns="http://www.opengis.net/kml/2.2")

    # Create <Document> as a child of <kml>
    document = Element_Tree.SubElement(kml, "Document")

    # Create a <Style> for the line (optional)
    style = Element_Tree.SubElement(document, "Style", id="lineStyle")
    line_style = Element_Tree.SubElement(style, "LineStyle")
    color = Element_Tree.SubElement(line_style, "color")
    color.text = "#ff0dde1e"  # color (AABBGGRR format)
    width = Element_Tree.SubElement(line_style, "width")
    width.text = "2"  # Line width

    # Iterate through the coordinates and create <Placemark> for each
    for i, (waypoint_name, (lat, lon, alt)) in enumerate(coordinates.items()):

        placemark = Element_Tree.SubElement(document, "Placemark")

        # Add <name> for each placemark (optional)
        name = Element_Tree.SubElement(placemark, "name")
        name.text = f"{waypoint_name}"

        # Set the scale of the icon to 0 to make it invisible
        icon_style = Element_Tree.SubElement(placemark, "Style")
        icon_style = Element_Tree.SubElement(icon_style, "IconStyle")
        scale = Element_Tree.SubElement(icon_style, "scale")
        scale.text = "0"  # Set scale to 0 to hide the icon

        # Set transparency (Alpha = 00 for full transparency)
        color = Element_Tree.SubElement(icon_style, "color")
        color.text = "00000000"  # Fully transparent color (AABBGGRR format)

        # Create <Point> element with <coordinates>
        point = Element_Tree.SubElement(placemark, "Point")
        coord = Element_Tree.SubElement(point, "coordinates")
        coord.text = f"{lon},{lat},0"  # KML uses "lon,lat,alt" format

    # Add a LineString connecting the points
    line_placemark = Element_Tree.SubElement(document, "Placemark")
    line_name = Element_Tree.SubElement(line_placemark, "name")
    line_name.text = "Path between Locations"
    style = Element_Tree.SubElement(line_placemark, "styleUrl")
    style.text = "#lineStyle"  # Use the previously defined line style

    # Create <LineString> and its coordinates
    line_string = Element_Tree.SubElement(line_placemark, "LineString")
    line_coordinates = Element_Tree.SubElement(line_string, "coordinates")
    # Join the coordinates into a single string: lon,lat,alt for each point
    line_coordinates.text = " ".join([f"{lon},{lat},{alt}" for (lat, lon, alt) in coordinates.values()])

    # Convert the tree structure to a string
    tree = Element_Tree.ElementTree(kml)

    # Create timestamp for the output file name
    current_time = datetime.now(timezone.utc)
    timestamp = " " + current_time.strftime("%H%M")
    file_name = output_label + timestamp + "Z.kml"

    # Write the KML to the output file
    with open(file_name, "wb") as f:
        tree.write(f, encoding="utf-8", xml_declaration=True)

    print("\n" + file_name + " created successfully.")
    print("If you are creating many .kml files, it is recommended to sort the file explorer by date/time to quickly "
          "identify the most recently generated file.\n")


# This check ensures that the script is being called from it's stored location, and not elsewhere.
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please send the details of this error to gordon.s.kiesling.mil@us.navy.mil")
        input("Press Enter to exit...")

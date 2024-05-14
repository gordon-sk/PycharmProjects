# This file & concept were created by LT Gordon "Yorkie" Kiesling, VUP-19 in May 2024.
#
# Go To Direct is a proposed MTD function for the MQ-4C Triton by Northrop Grumman.
# Please see the associated PowerPoint file "MTD GoToDirect Proposal.pptx" for further documentation
#
# GLOBAL RULES:
# This file was created in Pycharm Community Edition using Python 3.12.
# All calculations in this file are performed in radians, not degrees.
# For any reader and all pilots who are rusty on geometry, 360 degrees == 2*pi radians.
# All speeds are in meters per second.
# For the purposes of printing, everything in print commands is converted to knots or degrees.
#
# Libraries imported:
# MATH -- for calculations. Documentation on math: https://docs.python.org/3/library/math.html
# RANDOM -- for generating simulated parameters. In real life, the UA and/or AVP would provide everything that the
# random library provides here. Documentation on random: https://docs.python.org/3/library/random.html
# FOLIUM -- a random library of pretty mapping scripts I used to visualize results. In real life, the MTD will fill this
# role. Documentation on folio is available at https://python-visualization.github.io/folium/latest/index.html. Folio
# can be installed in python by executing "pip install folio" in a command shell or terminal.

import math
import random

try:
    import folium
except ModuleNotFoundError:
    print("\nNOTE: Folium is not installed. Proceeding with the rest of the program. All calculations will be " +
          "available to the end user, but graph products will not.")


def main():
    # GoToDirect is the master and controlling function for going to a waypoint directly, as described in the
    # PowerPoint file associated with this script.
    # INPUT: desired_waypoint is a 2-value library taking the form of [latitude, longitude]

    print("\n\nStart script!!\n")
    # BEGIN SIMULATION PARAMETERS -- defining values required for calculation.
    # IRL, many of these would be pulled from 1Hz and 1/5Hz data provided by the UA, and then obtained from the
    # Linux-run PFD software to the MTD.
    # A random heading anywhere on the compass:
    start_heading = math.radians((random.uniform(0, 360)))
    # defining boundaries for our initial, random lat/long:
    lat_position_floor = 16
    lat_position_ceiling = 16.5
    long_position_floor = 133
    long_position_ceiling = 133.5
    # position takes the form of [lat, long] and is stored in units of radians:
    position = [math.radians(random.uniform(lat_position_floor, lat_position_ceiling)),
                math.radians(random.uniform(long_position_floor, long_position_ceiling))]
    desired_waypoint = [math.radians(random.uniform(lat_position_floor, lat_position_ceiling)),
                        math.radians(random.uniform(long_position_floor, long_position_ceiling))]
    # wind takes the form of [direction (radians), speed (m/s)]:
    wind_speed_low = 12  # about 25 knots
    wind_speed_high = 39  # about 75 knots
    winds = [math.radians(random.uniform(0, 360)), (random.uniform(wind_speed_low, wind_speed_high))]
    # True Airspeed is variable, and can be between roughly 140 to 300 knots:
    TAS_low = 72  # m/s
    TAS_high = 155
    TAS = random.uniform(TAS_low, TAS_high)
    # turn rate -- randomly selected, but eventually can be calculated from NATOPS:
    turn_rate_low = 1.3  # deg/sec
    turn_rate_high = 3.1
    turn_rate = math.radians(random.uniform(turn_rate_low, turn_rate_high))  # converted to radians
    # misc -- unused for now:
    altitude = 46, 000
    weight = 30, 000

    # Some data validation printouts, for the author:
    print("The UA's current position is:\t\t\t\t\t\t\t" + c_string(position[0]) + "N, " + c_string(position[1]) + "E")
    print("The desired waypoint position is:\t\t\t\t\t\t" + c_string(desired_waypoint[0]) + "N, " + c_string(
        desired_waypoint[1]) + "E")
    print("Winds are currently:\t\t\t\t\t\t\t\t\t" + c_string(winds[0], 0) + " deg at " + knots(winds[1]))
    print("TAS at start of turn:\t\t\t\t\t\t\t\t\t" + knots(TAS))
    print("At it's current altitude, the UA will turn at:\t\t\t" + c_string(turn_rate, 1) + " deg/sec\n")

    # Calling supporting functions.
    new_heading = calculate_heading(position, desired_waypoint, TAS=TAS, winds=winds)
    turn_data = calculate_direction_of_turn(new_heading, start_heading)
    rough_heading_change = turn_data[0]
    direction_of_turn = turn_data[1]
    time_in_turn = rough_heading_change / turn_rate

    print("Turning for " + c_string(rough_heading_change, 0) + " degrees at " + str(
        round(math.degrees(turn_rate), 1)) + " deg/sec will take:\t\t" + str(round(time_in_turn, 1)) + " sec")

    # calculating the final, desired product -- an intermediate fix that will turn the UA towards the actual desired
    # waypoint.
    positions_list = calculate_intermediate_fix(time_in_turn, winds, TAS, start_heading, turn_rate, direction_of_turn)
    positions_list = convert_displacements(position, positions_list)
    intermediate_waypoint = positions_list[positions_list.__len__() - 1]
    print("\n\nIntermediate waypoint is:\t\t\t\t\t\t\t\t" + c_string(intermediate_waypoint[0], 4) + "N, " + c_string(
        intermediate_waypoint[1], 4) + "E")
    print("UA will be on a heading of:\t\t\t\t\t\t\t\t\t\t" + c_string(intermediate_waypoint[2], d=0))

    # Folium is the only library required to be installed -- and I really doubt NMCI will let us download and
    # install it via Pip. So, anticipating this script making it to that network, I've attempted to keep it
    # functional.
    try:
        plot(position, intermediate_waypoint, desired_waypoint, positions_list, lat_position_floor, long_position_floor,
             lat_position_ceiling, long_position_ceiling)
    except ModuleNotFoundError:
        print("Skipping Folium plotting portion of the script")

    print("End script")


def calculate_heading(current_pos, desired_waypoint, TAS=None, winds=None):
    # Source for initial course calculation: https://www.movable-type.co.uk/scripts/latlong.html
    # Source for wind correction: https://www.omnicalculator.com/physics/wind-correction-angle
    # INPUTS:
    # current_pos = [lat, long] of UA in radians at time of command execution
    # desired_waypoint = [lat, long] of a different waypoint
    # winds = [direction (radians), speed (knots)]. OPTIONAL
    # TAS = a random float representing UA True Airspeed at time of command execution. OPTIONAL
    # OUTPUT:
    # Returns heading, a float representing the azimuth/course/heading between two lat/long pairings. If winds are
    # included in the function arguments, a wind-corrected heading is returned.
    #
    # NOTE: This is an imperfect solution, as the actual new heading to fly will have changed due to the lateral
    # displacement of the UA that occurs during its turn. This error will be more pronounced during shorter-distance
    # Go To commands. For the purposes of this exercise I found the trigonometry too intimidating and I thus declare
    # myself satisfied with this mostly realized solution. My degree is in physics and my professors always told me
    # that you engineering types endorse this concept of "good enough" habitually ;)

    # splitting both lat/longs into components

    theta_1 = current_pos[0]
    lambda_1 = current_pos[1]
    theta_2 = desired_waypoint[0]
    lambda_2 = desired_waypoint[1]

    if winds is None:
        winds = [0, 0]
    if TAS is None:
        TAS = 0

    # solving for the course between the two sets of waypoints
    Y = math.sin(lambda_2 - lambda_1) * math.cos(theta_2)
    X = (math.cos(theta_1) * math.sin(theta_2)) - (
                math.sin(theta_1) * math.cos(theta_2) * math.cos(lambda_2 - lambda_1))
    azimuth = math.atan2(Y, X)
    azimuth = (azimuth + (2 * math.pi)) % (2 * math.pi)

    print("The azimuth (the course between two points) is:\t\t\t" + c_string(azimuth, 0) + " degrees TRUE")

    # calculate the crab angle as a result of winds
    delta_angle = (math.asin(winds[1] / TAS * math.sin((winds[0] - azimuth))))
    # application of crab angle, and an out-of-bounds sanity catch / normalization
    azimuth = (azimuth + delta_angle + (2 * math.pi)) % (2 * math.pi)
    print("The crab angle will be:\t\t\t\t\t\t\t\t\t" + str(round(math.degrees(delta_angle), 1)) + " degrees")
    print("The wind-corrected heading to fly is:\t\t\t\t\t" + c_string(azimuth, 0) + " degrees TRUE\n")
    return azimuth


def calculate_direction_of_turn(new_heading, current_heading):
    # This function is intended to calculate the direction of magnitude of a turn.
    # It accepts two inputs:
    # new_heading is the final rollout heading the AVP desires the UA to turn to.
    # current_heading is the heading of the UA at time of command execution.

    print("My current heading is:\t\t\t\t\t\t\t\t\t" + c_string(current_heading, 0) + " degrees")
    # For sanity checking purposes, I calculate turns in both directions -- effect on script runtime is negligible
    # out-of-bounds checking is done in the same line.
    right_turn_mag = (((2 * math.pi) - current_heading) + new_heading) % (2 * math.pi)
    left_turn_mag = (current_heading + ((2 * math.pi) - new_heading)) % (2 * math.pi)
    print("Right turn magnitude is:\t\t\t\t\t\t\t\t" + c_string(right_turn_mag, 0) + " degrees")
    print("Left turn magnitude is:\t\t\t\t\t\t\t\t\t" + c_string(left_turn_mag, 0) + " degrees")
    direction_of_turn = "LEFT"
    if right_turn_mag <= left_turn_mag:
        direction_of_turn = "RIGHT"
    print("Shorter direction of turn is:\t\t\t\t\t\t\t" + direction_of_turn)
    shorter_turn_magnitude = min(right_turn_mag, left_turn_mag)
    return [shorter_turn_magnitude, direction_of_turn]


def calculate_intermediate_fix(time_in_turn, winds, TAS, start_heading, turn_rate, direction_of_turn):
    # Imperfect solution bc it does not account for roll-in or roll-out of the turn
    # It also assumes a constant TAS, instead of a variable one as winds shift

    # latitudinal and longitudinal speed of wind
    wind_x = winds[1] * math.cos(winds[0])
    wind_y = winds[1] * math.sin(winds[0])
    # We define the instantaneous current position as 0m, 0m.
    pos_x = 0
    pos_y = 0
    # a new heading variable, which will be iterated over
    heading = start_heading
    # direction of turn: positive for a right turn, negative for a left turn
    direction_of_turn_variable = 1
    if direction_of_turn == "LEFT":
        direction_of_turn_variable = -1  # we have to go back Marty!
    # initial time, step size, and calculating how many times to iterate the loop
    t = 0
    step_size = .01
    number_of_iterations = int(time_in_turn / step_size)
    # saving the successive positions and headings for plotting later -- only for proof of concept
    positions_record = []

    # Performing the iterative calculation.
    for T in range(0, number_of_iterations):
        # First we split our TAS into latitudinal and longitudinal components.
        # Then, we add the respective components of wind, which are constant.
        V_x = TAS * math.cos(heading) + wind_x
        V_y = TAS * math.sin(heading) + wind_y
        # Next, we calculate a new position based on our old position, component velocity, and elapsed time.
        # distance = velocity (m/s) * t (seconds)
        pos_x = pos_x + (V_x * step_size)
        pos_y = pos_y + (V_y * step_size)

        # UA heading is updated based on rate of turn
        heading = heading + (turn_rate * step_size * direction_of_turn_variable)
        t += step_size
        positions_record.append([pos_x, pos_y, heading])
    if heading <= 0:
        heading += (2 * math.pi)  # all units in radians
    if heading >= (2 * math.pi):
        heading = heading % (2 * math.pi)
    print("\nAt completion of turn:")
    print("Elapsed time = " + str(round(t, 1)) + " seconds")
    print("UA position is " + str(round(pos_x / 1000, 1)) + " km lateral of starting position and " + str(
        round(pos_y / 1000, 1)) + " km longitudinal of starting position")
    print("UA heading at conclusion of turn is: " + c_string(heading, 0))

    return positions_record


def convert_displacements(original_position, positions_list):
    # We use distance (m) = radius (m) * omega (radians)
    # distance is our displacement, radius is that of the Earth, omega is what we solved for -- a tiny angle that
    # will give us theta_2 or lambda_2, IE a lat/long pair corresponding to some lateral and longitudinal displacement
    # from our origin point.
    # Solved, it takes the form of omega = radius / distance.
    # Because lines of latitude run parallel to the equator, we can simply use the Earth's radius for our calculations.
    # But lines of longitude run between the poles, and the Earth is variably smaller the closer you get to those poles.
    # Thus, we must calculate radius inside the loop for each longitude calculation.

    theta_1 = original_position[0]
    lambda_1 = original_position[1]
    earth_radius = 6371000  # in meters

    for a in range(0, positions_list.__len__()):
        x_displacement = positions_list[a][0]
        y_displacement = positions_list[a][1]
        theta_2 = theta_1 + (x_displacement / earth_radius)
        circumference = earth_radius * math.cos(theta_1)
        lambda_2 = lambda_1 + (y_displacement / circumference)
        positions_list[a][0] = theta_2
        positions_list[a][1] = lambda_2

    return positions_list


def plot(position, intermediate_waypoint, desired_waypoint, positions_list, lat_position_floor, long_position_floor, lat_position_ceiling,
         long_position_ceiling):
    # Finally -- plotting the finished product, and displaying calculated data. Don't dig too far into this, since
    # it won't be an option for the actual C# product
    # tiles is the URL by which Folium retrieves the desired map layer for display
    tiles = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
    # attribution is just good manners
    attribution = (
        'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community')
    filename = "UA_plot.html"
    # Initial map call
    map_file = folium.Map(
        location=(math.degrees(position[0]), math.degrees(position[1])),
        control_scale=True,
        zoom_start=8,
        tiles=tiles,
        attr=attribution
    )
    # adding a nice box that contains our random starting and ending position boundaries
    boundaries = folium.PolyLine(
        locations=[[lat_position_floor, long_position_floor],
                   [lat_position_ceiling, long_position_floor],
                   [lat_position_ceiling, long_position_ceiling],
                   [lat_position_floor, long_position_ceiling],
                   [lat_position_floor, long_position_floor]],
        color='red')
    plane_icon = folium.CustomIcon(
        'plane-up-solid.png',
        icon_anchor=(10, 10),
    )
    initial_pos_marker = folium.Marker(
        location=[math.degrees(position[0]), math.degrees(position[1])],
    )
    intermediate_pos_marker = folium.Marker(
        location=[math.degrees(intermediate_waypoint[0]), math.degrees(intermediate_waypoint[1])]
    )
    final_pos_marker = folium.Marker(
        location=[math.degrees(desired_waypoint[0]), math.degrees(desired_waypoint[1])]
    )

    coordinates = []
    for x in range(0, positions_list.__len__() - 1):
        coordinates.append(
            [math.degrees(positions_list[x][0]),
             math.degrees(positions_list[x][1])]
        )
    turn_line = folium.PolyLine(
        locations=coordinates,
        color='purple',
        weight=4,
        smooth_factor=0,
    )
    intial_course_line = folium.PolyLine(
        locations=(
            [math.degrees(position[0]), math.degrees(position[1])],
            [math.degrees(desired_waypoint[0]), math.degrees(desired_waypoint[1])]
        ),
        color='red',
        weight=1,
        smooth_factor=25,
    )
    final_course_line = folium.PolyLine(
        locations=(
            [math.degrees(intermediate_waypoint[0]), math.degrees(intermediate_waypoint[1])],
            [math.degrees(desired_waypoint[0]), math.degrees(desired_waypoint[1])]
        ),
        color='green',
        weight=1,
        smooth_factor=25,
    )

    boundaries.add_to(map_file)
    turn_line.add_to(map_file)
    intial_course_line.add_to(map_file)
    final_course_line.add_to(map_file)
    initial_pos_marker.add_to(map_file)
    intermediate_pos_marker.add_to(map_file)
    final_pos_marker.add_to(map_file)

    map_file.save(filename)
    open(filename)


def c_string(n, d=2):
    # This function accepts a float as input, with an option parameter for desired decimal place to round the float.
    # It is intended to accept a measurement in radians, convert it to degrees, and convert that to a string for
    # printing. Decimal place can be optionally specified and defaults to the hundredths place.

    n_deg = (math.degrees(n) + 360) % 360

    final = str(round(n_deg, d))
    # conditional formatting to make the printout pilot-friendly
    if d == 0:
        final = str(int(round(n_deg, d)))
        if n_deg < 10:
            final = "00" + str(int(round(n_deg, 0)))
        else:
            if n_deg < 100:
                final = "0" + str(int(round(n_deg, 0)))
    return final


def knots(n):
    # accepts n (m/s), returns knots in a printable form
    return str(int(round(n * 1.9438444924, 0))) + " knots"


main()

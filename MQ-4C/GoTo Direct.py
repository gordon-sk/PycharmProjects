# This file & concept were created by LT Gordon "Yorkie" Kiesling, VUP-19 in May 2024.
#
# Go To Direct is a proposed MTD function for the MQ-4C Triton by Northrop Grumman.
# Please see the associated PowerPoint file "MTD GoToDirect Proposal.pptx" for further background.
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
# can be installed in python by executing "pip install folio" in a command shell or terminal. If unable, the script
# will function without it.

import math
import random
import os

try:
    import folium
except ModuleNotFoundError:
    print("\nNOTE: Folium is not installed. Proceeding with the rest of the program. All calculation results will be " +
          "available to you, the visual product will not.")


def main():
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
    start_position = [math.radians(random.uniform(lat_position_floor, lat_position_ceiling)),
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
    # altitude = 46,000
    # weight = 30,000

    # Some data printouts:
    print("The UA's current position is:\t\t\t\t\t\t\t" + c_string(start_position[0]) + "N, " +
          c_string(start_position[1]) + "E")
    print("The desired waypoint position is:\t\t\t\t\t\t" + c_string(desired_waypoint[0]) + "N, " + c_string(
        desired_waypoint[1]) + "E")
    print("Winds are currently:\t\t\t\t\t\t\t\t\t" + c_string(winds[0], 0) + " deg at " + knots(winds[1]))
    print("TAS at start of turn:\t\t\t\t\t\t\t\t\t" + knots(TAS))
    print("Heading at start of turn:\t\t\t\t\t\t\t\t" + c_string(start_heading, 0) + " degrees TRUE")
    print("At it's current altitude, the UA will turn at:\t\t\t" + c_string(turn_rate, 1) + " deg/sec\n")

    # Calling supporting functions.
    rough_new_heading = calculate_azimuth(start_position, desired_waypoint, TAS=TAS, winds=winds)
    print("At start position, wind-corrected azimuth is:\t\t\t" + c_string(rough_new_heading, d=0) + " degrees TRUE")

    short_turn_magnitude, direction_of_turn = calculate_direction_of_turn(rough_new_heading, start_heading)
    print("The turn should take roughly:\t\t\t\t\t\t\t" + str(round(short_turn_magnitude/turn_rate, 1)) + " sec")

    # calculating the final, desired product -- an intermediate fix that will turn the UA towards the actual desired
    # waypoint. Actually, we get a list of lat/longs, because we want to plot them -- but IRL this would simply return
    # a single lat/long for a Goto command execution.
    positions_list = calculate_intermediate_fix(start_position, desired_waypoint, winds, TAS, start_heading,
                                                rough_new_heading, turn_rate, direction_of_turn)

    # Pulling and printing the final lat/long from the calculation, which is what the MTD would actually issue a
    # GoTo command for. The follow-on waypoint would be the desired_waypoint.
    intermediate_waypoint = positions_list[positions_list.__len__() - 1]
    print("\nintermediate waypoint is:\t\t\t\t\t\t\t\t" + c_string(intermediate_waypoint[0], 4) + "N, " + c_string(
        intermediate_waypoint[1], 4) + "E\n")

    # Folium is the only library that doesn't come default with Python -- and I really doubt NMCI will let us download
    # and install it via Pip. So, anticipating this making it to that network, I've attempted to keep the script
    # functional without Folium present.
    try:
        plot(start_position, intermediate_waypoint, desired_waypoint, positions_list, lat_position_floor,
             long_position_floor, lat_position_ceiling, long_position_ceiling)
    except NameError:
        print("Skipping Folium plotting portion of the script -- visual product unavailable")

    print("End script")


def calculate_azimuth(start_pos, desired_waypoint, TAS=None, winds=None):
    # Source for initial course calculation: https://www.movable-type.co.uk/scripts/latlong.html
    # Source for wind correction: https://www.omnicalculator.com/physics/wind-correction-angle
    # This function calculates the azimuth between two latitude/longitude pairs. Wind correction is optional, and if
    # elected requires both TAS and wind data.
    # INPUTS:
    # current_pos = [lat, long] of UA in radians at time of command execution
    # desired_waypoint = [lat, long] of a different waypoint
    # winds = [direction (radians), speed (knots)]. OPTIONAL
    # TAS = a float representing UA True Airspeed at time of command execution. OPTIONAL
    # OUTPUT:
    # Returns heading, a float representing the azimuth/course/heading between two lat/long pairings. If winds are
    # included in the function arguments, a wind-corrected heading is returned.

    # splitting both lat/longs into components.
    theta_1 = start_pos[0]
    lambda_1 = start_pos[1]
    theta_2 = desired_waypoint[0]
    lambda_2 = desired_waypoint[1]

    # solving for the course between the two sets of waypoints.
    Y = math.sin(lambda_2 - lambda_1) * math.cos(theta_2)
    X = (math.cos(theta_1) * math.sin(theta_2)) - (
                math.sin(theta_1) * math.cos(theta_2) * math.cos(lambda_2 - lambda_1))
    azimuth = math.atan2(Y, X)
    azimuth = (azimuth + (2 * math.pi)) % (2 * math.pi)

    if winds is not None and TAS is not None:
        # calculate the crab angle as a result of winds.
        delta_angle = (math.asin(winds[1] / TAS * math.sin((winds[0] - azimuth))))
        # application of crab angle, and an out-of-bounds sanity catch / normalization.
        azimuth = (azimuth + delta_angle + (2 * math.pi)) % (2 * math.pi)
    else:
        print("Wind correction requires both TAS and winds. Wind correction skipped because one of those two missing.")

    return azimuth


def calculate_direction_of_turn(new_heading, current_heading):
    # This function is intended to calculate the direction of magnitude of a turn.
    # INPUTS:
    # new_heading is the final rollout heading the AVP desires the UA to turn to.
    # current_heading is the heading of the UA at time of command execution.
    # OUTPUTS:
    # shorter_turn_magnitude: a float, in radians, of the estimate of the size of turn required to accomplish the
    # passed heading change.
    # direction_of_turn: A string, either "LEFT" or "RIGHT"

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
    return shorter_turn_magnitude, direction_of_turn


def calculate_intermediate_fix(start_pos, desired_waypoint, winds, TAS, start_heading, rough_new_heading,
                               turn_rate, direction_of_turn):
    # The meat and potatoes of our work here!
    # This function takes our starting position, desired waypoint, air data, and turn estimations and iterates by step
    # size until we expect the turn to be completed. At each step, it determines a new position and heading for the UA.
    # INPUTS:
    # start_pos: The UA's position at time of command execution, formatted as [lat, long].
    # desired_waypoint: The [lat, long] of the AVP specified desired waypoint.
    # winds: wind data at time of command execution, formatted as [wind direction, wind speed] in [radians, m/s].
    # TAS: True airspeed at time of command execution in m/s.
    # start_heading: UA heading in radians at time of command execution.
    # rough_new_heading: an estimate of rollout heading at time of command execution.
    # turn_rate: instantaneous expected turn rate of the UA at time of command execution, in radians/sec.
    # direction_of_turn: either "LEFT" or "RIGHT" as previously determined.
    # OUTPUT: positions_record, a library of libraries where the latter takes the form of [lat, long, heading] at each
    # iteration point in our while loop.

    # initializing time variable.
    t = 0  # seconds
    # choosing step size to iterate with.
    step_size = .01  # seconds

    # latitudinal and longitudinal speed of wind.
    wind_x = winds[1] * math.cos(winds[0])
    wind_y = winds[1] * math.sin(winds[0])
    # We define the instantaneous current position as 0m, 0m.
    x_displacement = 0
    y_displacement = 0
    # a new heading variable, which will be updated as the UA turns. We also initialize a variable to store the current
    # "planned" rollout heading, and the remaining heading change to turn.
    current_heading = start_heading
    current_rollout_heading = rough_new_heading
    outstanding_heading_change = math.fabs(current_heading - current_rollout_heading)
    # direction of turn: positive for a right turn, negative for a left turn
    direction_of_turn_variable = 1
    if direction_of_turn == "LEFT":
        direction_of_turn_variable = -1  # we have to go back Marty!

    # saving the successive positions and headings for plotting later.
    positions_record = []

    # Performing the iterative calculation. We iterate until the turn remaining is less than 1 degree, which is more
    # than good enough for a step size of .01.
    iteration_counter = 0
    while outstanding_heading_change > math.radians(1.0):
        # First we split our TAS into latitudinal and longitudinal components.
        # Then, we add the respective components of wind, which we consider constant.
        V_x = TAS * math.cos(current_heading) + wind_x
        V_y = TAS * math.sin(current_heading) + wind_y
        # Next, we calculate a new position based on our old position, component velocity, and elapsed time.
        # distance = velocity (meters / second) * t (seconds)
        x_displacement = x_displacement + (V_x * step_size)
        y_displacement = y_displacement + (V_y * step_size)
        # Calling our conversion function to change the displacements to lat longs, relative to UA start position.
        updated_lat, updated_long = convert_displacement(start_pos, x_displacement, y_displacement)

        # We calculate a new rollout heading, which will have changed due to the position of the UA updating.
        new_rollout_heading = calculate_azimuth(start_pos=[updated_lat, updated_long],
                                                desired_waypoint=desired_waypoint,
                                                TAS=TAS,
                                                winds=winds)
        # With a new rollout heading comes a new "amount of turn" remaining. With sanity checks.
        outstanding_heading_change = math.fabs(current_heading-new_rollout_heading)
        if outstanding_heading_change < 0:
            print("Negative heading change remaining -- this shouldn't be possible. You may have made step size too"
                  "large. The results of... whatever you did, will appear below.")
            break
        if outstanding_heading_change >= (2 * math.pi):
            outstanding_heading_change = outstanding_heading_change % (2 * math.pi)

        # Current UA heading is updated based on rate of turn. With sanity checks.
        current_heading = current_heading + (turn_rate * step_size * direction_of_turn_variable)
        if current_heading <= 0:
            current_heading += (2 * math.pi)
        if current_heading >= (2 * math.pi):
            current_heading = current_heading % (2 * math.pi)

        # Recording time elapsed, iterations, and new UA position. Loop continues until the turn is complete.
        t += step_size
        iteration_counter += 1
        positions_record.append([updated_lat, updated_long, current_heading])

    # Some fun data printout for the end user... IE my sanity.
    print("\nAt completion of turn:")
    print("Elapsed time:\t\t\t\t\t\t\t\t\t\t\t" + str(round(t, 1)) + " seconds")
    print("Elapsed iterations:\t\t\t\t\t\t\t\t\t\t" + str(iteration_counter))
    print("UA lateral displacement is\t\t\t\t\t\t\t\t" + str(round(x_displacement / 1000, 1)) + " km")
    print("UA longitudinal displacement is\t\t\t\t\t\t\t" + str(round(y_displacement / 1000, 1)) + " km")
    print("UA heading at conclusion of turn is:\t\t\t\t\t" + c_string(current_heading, 0) + " degrees TRUE")

    return positions_record


def convert_displacement(start_pos, x_displacement, y_displacement):
    # We use distance (m) = radius (m) * omega (radians)
    # distance is our displacement, radius is that of the Earth, omega is what we solved for -- a tiny angle that
    # will give us theta_2 or lambda_2, IE a lat/long pair corresponding to some lateral and longitudinal displacement
    # from our origin point.
    # Solved, it takes the form of omega = radius / distance.
    # Because lines of latitude run parallel to the equator, we can simply use the Earth's radius for our calculations.
    # But lines of longitude run between the poles, and the Earth is variably smaller the closer you get to those poles.
    # Thus, we must calculate radius inside the loop for each longitude calculation. This variable is called
    # circumference below.
    # INPUT:
    # start_pos: a [lat, long] pair intended to be the start position of the UA at time of command execution.
    # x_displacement: A float representing the displacement of the UA in meters laterally from start position.
    # y_displacement: Same same, but longitudinal.
    # OUTPUT:
    # lat, long corresponding to the displacement from start position in meters.

    theta_1 = start_pos[0]
    lambda_1 = start_pos[1]
    earth_radius = 6371000  # in meters

    theta_2 = theta_1 + (x_displacement / earth_radius)
    circumference = earth_radius * math.cos(theta_1)
    lambda_2 = lambda_1 + (y_displacement / circumference)

    return theta_2, lambda_2


def plot(position, intermediate_waypoint, desired_waypoint, positions_list, lat_position_floor, long_position_floor,
         lat_position_ceiling, long_position_ceiling):
    # Finally -- plotting the finished product, and displaying calculated data. Don't dig too far into this, since
    # it won't be an option for the actual C# product
    # tiles is the URL by which Folium retrieves the desired map layer for display
    tiles = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
    # attribution is just good manners
    attribution = (
        'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, '
        'UPR-EGP, and the GIS User Community')
    filename = "UA_plot.html"
    # Initial map call
    map_file = folium.Map(
        location=(math.degrees(position[0]), math.degrees(position[1])),
        control_scale=True,
        zoom_start=12,
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

    coordinates = []
    plane_locations = []
    plane_headings = []
    for x in range(0, positions_list.__len__() - 1):
        coordinates.append(
            [math.degrees(positions_list[x][0]),
             math.degrees(positions_list[x][1])]
        )
        # Every 300 points, we save coordinates so that we can draw a little plane and heading line
        if x % 300 == 0:
            plane_locations.append(
                [math.degrees(positions_list[x][0]),
                 math.degrees(positions_list[x][1])]
            )
            heading = (math.degrees(positions_list[x][2]) + 360) % 360
            plane_headings.append(heading)

    plane1_icon = folium.Icon(
        icon="plane",
        color="pink",
        angle=int((math.degrees(positions_list[0][2]) + 360) % 360)
    )
    plane2_icon = folium.Icon(
        icon="plane",
        color="green",
        angle=int((math.degrees(positions_list[positions_list.__len__()-1][2]) + 360) % 360)
    )
    desired_waypoint_icon = folium.Icon(
        icon="location-crosshairs"
    )

    initial_pos_marker = folium.Marker(
        location=[math.degrees(position[0]), math.degrees(position[1])],
        icon=plane1_icon,
    )
    intermediate_pos_marker = folium.Marker(
        location=[math.degrees(intermediate_waypoint[0]), math.degrees(intermediate_waypoint[1])],
        icon=plane2_icon
    )
    final_pos_marker = folium.Marker(
        location=[math.degrees(desired_waypoint[0]), math.degrees(desired_waypoint[1])],
        icon=desired_waypoint_icon
    )

    turn_line = folium.PolyLine(
        locations=coordinates,
        color='purple',
        weight=4,
        smooth_factor=0,
        popup="Expected flight path in turn"
    )
    initial_course_line = folium.PolyLine(
        locations=(
            [math.degrees(position[0]), math.degrees(position[1])],
            [math.degrees(desired_waypoint[0]), math.degrees(desired_waypoint[1])]
        ),
        color='red',
        weight=1,
        smooth_factor=25,
        popup="Initial IMMC calculated course"
    )
    final_course_line = folium.PolyLine(
        locations=(
            [math.degrees(intermediate_waypoint[0]), math.degrees(intermediate_waypoint[1])],
            [math.degrees(desired_waypoint[0]), math.degrees(desired_waypoint[1])]
        ),
        color='green',
        weight=1,
        smooth_factor=25,
        popup="Final, flown course"
    )

    boundaries.add_to(map_file)
    turn_line.add_to(map_file)
    initial_course_line.add_to(map_file)
    final_course_line.add_to(map_file)
    initial_pos_marker.add_to(map_file)
    intermediate_pos_marker.add_to(map_file)
    final_pos_marker.add_to(map_file)

    map_file.save(filename)
    os.startfile(filename)


def c_string(n, d=2):
    # This function accepts a float as input, with an option parameter for desired decimal place to round the float.
    # It is intended to accept a measurement in radians, convert it to degrees, and convert that to a string for
    # printing. Decimal place can be optionally specified and defaults to the hundredths place.
    # This is just to make stuff print friendly.

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
    # for ease of reading data printouts
    return str(int(round(n * 1.9438444924, 0))) + " knots"


main()

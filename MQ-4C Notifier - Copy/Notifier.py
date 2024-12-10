import time
from datetime import datetime, timezone, timedelta
import ctypes
import os
import sys
import math
import threading
import winsound

RED = "\033[31m"
YELLOW = "\033[33m"
BLACK_BG = "\033[40m"  # Black background
BRIGHT_GREEN_TEXT = "\033[92m"  # Bright green text
RESET = "\033[0m"  # Resets color back to default
os.system('color 0A')
hwnd = ctypes.windll.kernel32.GetConsoleWindow()
ack = False


def main():
    print(BLACK_BG + BRIGHT_GREEN_TEXT + "Welcome to the VUP-19 reminder tool.")
    print("Making sure your 2P is SOP compliant since 2024.")
    print("You can use this tool to remind you to take detailed statuses at predetermined time intervals.\n")
    start_time = get_first_reminder_time()
    winsound.Beep(659, 200)
    interval_seconds = get_interval()
    winsound.Beep(659, 200)
    label = get_label()
    winsound.Beep(659, 200)

    print("All set. You can minimize this window or push it to the background. It will pop up into the foreground "
          "automatically when needed.")
    winsound.Beep(880, 200)

    stop_event = threading.Event()

    while True:
        now = datetime.now(timezone.utc)
        reminder_time = start_time + timedelta(seconds=interval_seconds)
        if now > reminder_time:
            print()
            stop_event.clear()
            popup_thread = threading.Thread(target=reminder_popup, args=(label, stop_event))
            popup_thread.start()
            popup_ack(stop_event)
            start_time += timedelta(seconds=interval_seconds)
        else:
            seconds_remaining = (reminder_time - now).seconds
            if seconds_remaining >= 60:
                print("\rTime to reminder: " + YELLOW + str(int(math.floor(seconds_remaining/60))).zfill(2)
                      + BRIGHT_GREEN_TEXT + " minutes", end="", flush=True)
            else:
                print("\rTime to next reminder: " + YELLOW + str(seconds_remaining).zfill(2) + BRIGHT_GREEN_TEXT
                      + " seconds", end="", flush=True)
        time.sleep(1)


def get_first_reminder_time():
    while True:
        print("The current time is: " + RED + datetime.now(timezone.utc).strftime(
            "%H:%M") + "Z" + BLACK_BG + BRIGHT_GREEN_TEXT + ".")
        reminder_time_str = input("When would you like your first reminder? You can type " + RED + """now""" +
                                  BLACK_BG + BRIGHT_GREEN_TEXT + ", or enter a time in " + RED + "HHMM" +
                                  BLACK_BG + BRIGHT_GREEN_TEXT + " format: ")
        reminder_time_str = reminder_time_str.strip().lower()

        try:
            now = datetime.now(timezone.utc)
            if reminder_time_str == "now":
                print("First reminder will appear " + YELLOW + "ASAP" + BRIGHT_GREEN_TEXT + BLACK_BG + ".\n")
                return now
            else:
                if len(reminder_time_str) != 4:
                    raise ValueError("Input must be in HHMM format")

                hour = int(reminder_time_str[0:2])
                minute = int(reminder_time_str[2:4])

                if not (0 <= hour <= 23):
                    raise ValueError("Hour must be between 00 and 23")
                if not (0 <= minute <= 59):
                    raise ValueError("Minute must be between 00 and 59")


                else:
                    reminder_time = now.replace(hour=hour, minute=minute, second=0)
                    if reminder_time < now:
                        print("First reminder will appear tomorrow (ZULU) at: " + YELLOW +
                              f"{reminder_time.strftime('%H:%M')}Z\n" + BLACK_BG + BRIGHT_GREEN_TEXT)
                        return reminder_time + timedelta(days=1)
                    else:
                        print(f"First reminder will appear at: " + YELLOW + f"{reminder_time.strftime('%H:%M')}Z\n"
                              + BLACK_BG + BRIGHT_GREEN_TEXT)
                        return reminder_time

        except ValueError as e:
            print(RED + "ERROR " + BLACK_BG + BRIGHT_GREEN_TEXT + f"Incorrect formatting detected: {e}. "
                                                                  "Please try again.\n")


def get_interval():
    while True:
        print("What time interval would you like between reminders?")
        print("You can type " + RED + """ice""" + BLACK_BG + BRIGHT_GREEN_TEXT + " for 60 seconds, " + RED + """SOP"""
              + BLACK_BG + BRIGHT_GREEN_TEXT + " for 1 hour, or a custom number between " + RED + "1 - 1200"
              + BLACK_BG + BRIGHT_GREEN_TEXT + " for minutes, 1200 minutes being 20 hours. Useful for malfunctions "
              "where you want to keep an eye on a system, like an ECS 43.")
        reminder_time_interval = input("Interval: ")
        reminder_time_interval = reminder_time_interval.lower().strip()

        if reminder_time_interval == "ice":
            print("Ice mode. Reminder interval is set to " + YELLOW + "60 seconds" + BLACK_BG + BRIGHT_GREEN_TEXT
                  + ".\n")
            return 60
        elif reminder_time_interval == "sop":
            print("SOP mode. Reminder interval is set to " + YELLOW + "1 hour" + BLACK_BG + BRIGHT_GREEN_TEXT
                  + ".\n")
            return 60 * 60
        else:
            try:
                reminder_time_interval = int(reminder_time_interval)
                if not (1 <= reminder_time_interval <= 1200):
                    raise ValueError("Hour must be between 1 and 1200")
                print("Custom reminder interval set to " + YELLOW + str(reminder_time_interval) + " minute(s)" + BLACK_BG + BRIGHT_GREEN_TEXT + ".\n")
                return reminder_time_interval * 60

            except ValueError as e:
                print(RED + "ERROR: " + BLACK_BG + BRIGHT_GREEN_TEXT + f"Invalid interval format: {e}. "
                                                                       "Please try again.\n")


def get_label():
    print("Optionally, specify something to label the reminder with. For example, ""ice det"", ""ECS"". "
          "Hit enter to leave blank.")
    label = input("System: ")
    if label == "":
        print("Reminder label blank.\n")
        return ""
    else:
        print("Reminder label set to " + YELLOW + label + BLACK_BG + BRIGHT_GREEN_TEXT + ".\n")
        return label


def reminder_popup(system, stop_event):

    spin_chars = ['|', '/', '|', '\\']
    idx = 0

    winsound.Beep(660, 800)
    sys.stdout.write("\n\r\t\t\tEXECUTE detailed status: " + system + "\n")
    while not stop_event.is_set():

        if hwnd != 0:
            # Show the window
            ctypes.windll.user32.ShowWindow(hwnd, 9)

            # Bring the window to the foreground
            ctypes.windll.user32.SetForegroundWindow(hwnd)

        # print a bunch
        sys.stdout.write(YELLOW + "\r\t\t\tPress ENTER twice to reset" + BLACK_BG + BRIGHT_GREEN_TEXT)
        # Print the spinner character and overwrite the line
        sys.stdout.write("\r\t\t" + spin_chars[idx] + "\t\t\t\t\t" + spin_chars[idx])
        sys.stdout.flush()
        idx = (idx + 1) % len(spin_chars)

        time.sleep(0.1)

    print()


def popup_ack(stop_event):
    input()
    stop_event.set()
    winsound.Beep(880, 400)

main()

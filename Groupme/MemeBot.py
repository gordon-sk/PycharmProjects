access_token = 'Qvk2qpeSBvWTChmwGy6kH2m8RXOUnOL7sxrebcA2'

from groupy import Group, Bot
import time
import os

while True:
    try:
        groups = Group.list()
        for x in groups:
            if x.group_id == '32600517': # beta pi
                group = x
                break
            else:
                group = 0

        if group is not 0:
            for x in group.members():
                if x.id == '256865122': # ary
                    target = x
                    group.remove(target)
                    os.system("osascript -e \'display notification \"Ary Removed\" with title \"Nice\"")

    except:
        pass
    time.sleep(2)
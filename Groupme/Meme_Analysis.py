from __future__ import division
import matplotlib.pyplot as plt

# index 0 -- name of groupme
# index 1 -- date, format YEAR-MONTH-DAY TIME(MILITARY)
# index 2 -- name of poster
# index 3 -- message content, an empty string if just an image
# index 4 -- number of likes
filename = "morememes.csv"

def main():


    data = open_file(filename)
    print "Data read..."
    print data.__len__()
    post_count_dict = post_counting(data)
    likes_count_dict = like_counting(data)
    avg_likes_dict = avg_likes(post_count_dict, likes_count_dict)
    graphing(post_count_dict, likes_count_dict, avg_likes_dict)
    writing(post_count_dict, likes_count_dict, avg_likes_dict)


def post_counting(data):
    return_list = {}
    for x in data:
        return_list[x[3]] = 0
    for x in data:
        return_list[x[3]] = return_list[x[3]] + 1
    return return_list


def like_counting(data):
    return_list = {}
    for x in data:
        return_list[x[3]] = 0
    for x in data:
        return_list[x[3]] = return_list[x[3]] + int(x[5])
    return return_list


def avg_likes(posts, likes):
    return_list = {}
    for name in posts:
        try:
            return_list[name] = round(likes[name] / posts[name], 2)
        except ZeroDivisionError:
            return_list[name] = 0
    return return_list


def graphing(posts, likes, avgs):
    x_posts, y_likes = [], []
    for name in posts:
        x_posts.append(posts[name])
        y_likes.append(avgs[name])


    plt.plot(x_posts, y_likes, 'go')
    plt.savefig("graphed_data.png")


def writing(posts, likes, avgs):
    newfile = open(filename[:-4] + "_writeup.txt", 'w')
    newfile.write("Meme Analysis\n\n")
    newdict = {}
    for name in avgs:
        newdict[avgs[name]] = name
    print newdict

    for avg in sorted(newdict.keys(), reverse=True):
        name = newdict[avg]
        if posts[name] > 10:
            newfile.write("\nName: " + name)
            newfile.write("\nPost Count: " + str(posts[name]))
            newfile.write("\nTotal Likes: " + str(likes[name]))
            newfile.write("\nAvg likes/post: " + str(avgs[name]) + "\n\n")

    newfile.close()



def open_file(filename):
    data_file = open(filename)
    data = data_file.read().splitlines()
    return_list = []
    date_restriction = True
    if not date_restriction:
        for msg in data:
            sub = msg.split(',')
            if sub.__len__() == 6 and sub[4] == "":
                return_list.append(sub)
    else:
        for msg in data:
            sub = msg.split(",")
            if sub.__len__() == 6 and sub[4] == "":
                date = sub[2]
                date = date.split("-")
                if date[0] == "2017" and int(date[1]) >= 2:
                    return_list.append(sub)
    return return_list

main()
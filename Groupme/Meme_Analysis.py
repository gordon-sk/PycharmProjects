import matplotlib.pyplot as mp

# index 0 -- name of groupme
# index 1 -- date, format YEAR-MONTH-DAY TIME(MILITARY)
# index 2 -- name of poster
# index 3 -- message content, an empty string if just an image
# index 4 -- number of likes

def main():
    filename = "Memes.csv"
    data = open_file(filename)
    posts_per_person = counting_posts(data)
    likes_per_person = counting_likes(data)
    #removed = removal(data)
    likes_per_post_versus_number_of_posts(posts_per_person, likes_per_person)
    write_output(posts_per_person, likes_per_person, filename)


def likes_per_post_versus_number_of_posts(post_counts, like_counts):
    names, number_of_posts, likes_per_post_avg = [], [], []
    for x in post_counts:
        if post_counts[x] > 25 and x != "GroupMe":
            names.append(x)
            number_of_posts.append(post_counts[x])
            try:
                likes_per_post_avg.append(round((like_counts[x]*1.0)/post_counts[x], 2))
            except KeyError:
                pass
            except ZeroDivisionError:
                likes_per_post_avg.append(0)
    print number_of_posts.__len__()
    print likes_per_post_avg.__len__()
    print names.__len__()

    mp.plot(number_of_posts, likes_per_post_avg, 'ro')
    for x in range(names.__len__()-1):
        try:
            mp.annotate(names[x], xy=(number_of_posts[x], likes_per_post_avg[x]))
        except ValueError:
            pass
    mp.ylabel("avg likes per post")
    mp.xlabel("total number of posts")
    mp.show()


def open_file(filename):
    data_file = open(filename)
    readable = data_file.read().splitlines()
    data_file.close()    # One of my profs once told me if we didn't include this, it was sloppy code, so here it is
    memes, retard_data = [], []           # Our master list of data
    for x in range(readable.__len__()):     # Here we separate post data
        sub = readable[x].split(',')        # splitting space-separated data
        if sub.__len__() != 5:
            retard_data.append(sub)
        else:
            memes.append(sub)                   # Adding the sublist to our master list
    return memes


def counting_posts(data):
    names = dict()
    for x in data:
        if x.__len__() >= 3:
            if names.has_key(str(x[2])):
                names[str(x[2])] += 1
            else:
                names[str(x[2])] = 0
    return names


def counting_likes(data):
    names = dict()
    for x in data:
        if x.__len__() == 5:
            if names.has_key(str(x[2])):
                try:
                    names[str(x[2])] += int(x[4])
                except ValueError:
                    pass
            else:
                names[str(x[2])] = 0
    return names


def removal(data):
    removals = []
    for x in data:
        if x.__len__() >= 5:
            if x[3].__contains__("removed"):
                removals.append(x[3].split(' '))
    for x in removals:
        if x.__contains__('from') or x.__contains__('the') or x.__contains__('group'):
            x.remove("from")
            x.remove("the")
            x.remove("group")
    print removals
    return removals


def write_output(posts_per_person, likes_per_person, data_name):
    final_data = []
    for x in posts_per_person:
        try:
            posts = posts_per_person[x]
            likes_total = likes_per_person[x]
            try:
                avg_likes = round(1.0*likes_total/posts, 2)
            except ZeroDivisionError:
                avg_likes = 0
            final_data.append([x, posts, likes_total, avg_likes])
        except KeyError:
            pass

    avg_high = 0
    avg_low = 100000
    avg_guy = ''
    avg_shit = ''
    for x in final_data:
        if x[3] > avg_high:
            avg_high = x[3]
            avg_guy = x[0]
        elif 0 < x[3] <= avg_low and x[1] > 5 and x[0] != 'GroupMe':
            avg_low = x[3]
            avg_shit = x[0]

    output = open(data_name + "analysis", "w")
    output.write(avg_guy + " has the best average, with a value of " + str(avg_high) + " likes/post\n")
    output.write(avg_shit + str(" has the worst average, at " + str(avg_low) + " likes per post\n\n"))
    output.write("Name, Number of Posts, total likes, avg posts/like\n")
    for x in final_data:
        output.write(str(x) + "\n")
    output.close()


main()
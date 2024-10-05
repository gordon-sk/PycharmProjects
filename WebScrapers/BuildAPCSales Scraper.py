import praw
import time
from email.mime.text import MIMEText
import smtplib
master_seen_list = []


def scraper():
    reddit = praw.Reddit(client_id='rP2oAQNBpj5uOg',
                         client_secret="ZPHVbBZqAs-GMrUtNNRlm8gVU2A", password='!bOMe6u&vFP%Tw8z',
                         user_agent='USERAGENT', username='anon_scraper')


    try:
        subreddit = reddit.subreddit('buildapcsales')
        strings = ['Cooler', 'cooler', 'COOLER', '9600', '660p', '9700', 'z390', 'Z390']
        hits = []
        x=0
        for submission in subreddit.new():
            print(x)
            for text in strings:
                if text in submission.title.lower():
                    if master_seen_list.__contains__(submission):
                        break
                    else:
                        hits.append(submission)
                        master_seen_list.append(submission)
            x += 1
        if len(hits) > 0:
            email(hits)
            print("sending hits")
        else:
            print('None found')
    except Exception:
        print(Exception)
        pass


def email(hits):
    # Email code
    sender = 'burner.checker@gmail.com'
    receiver = 'G.s.kiesling@gmail.com'
    body = ''
    if hits.__len__() == 1:
        body += 'There is a new post on /r/buildapcsales regarding a computer part.\n\n'
    else:
        body += 'There are new posts on r/buildapcsales regarding computer parts\n\n'
    for hit in hits:
        body += str(hit.title) + ':     \nhttps://old.reddit.com/' + str(hit.permalink)
        body += '\n\n'

    # Structuring the email
    msg = MIMEText(body)
    msg['From'] = sender
    msg['To'] = receiver
    msg['Subject'] = "New Post on /r/buildapcsales"

    # Access sender email and send email
    mail = smtplib.SMTP('smtp.gmail.com', 587)
    mail.ehlo()
    mail.starttls()
    mail.login('burner.checker@gmail.com', 'Rhsk24fn')
    message = msg.as_string()
    mail.sendmail(sender, receiver, message)
    mail.close()
    print('hits sent')


if __name__ == '__main__':
    while True:
        scraper()
        time.sleep(5*60)

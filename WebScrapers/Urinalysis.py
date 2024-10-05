import time
from email.mime.text import MIMEText
import smtplib



def scraper():




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
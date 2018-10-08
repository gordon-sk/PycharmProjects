from bs4 import BeautifulSoup as bs
import time
import urllib.request
from email.mime import multipart
from email.mime.text import MIMEText
import smtplib

check, kickoff = False, None
url = 'http://www.vucommodores.com/sports/m-footbl/sched/vand-m-footbl-sched.html'
checklist = ['      vs. Western Kentucky', '      Nashville, Tenn.', '      11/04/17', '      TBA', '', '      ']

while True:
    response = urllib.request.urlopen(url)
    webContent = response.read()
    webContent = str(webContent)
    f = open('vandyfootball.html', 'w')
    f.write(webContent)
    f.close

    file = open('vandyfootball.html')
    soup = bs(file, 'html.parser')
    lines = soup.text.split('\\n')
    for x in lines:
        if x == '      11/04/17':
            check = True
        if x == 'Audio':
            check = False
        if check:
            if checklist.__contains__(x):
                pass
            else:
                print(x)
                kickoff = x

    if kickoff is not None:
        # Email code
        sender = 'aa.http.check@gmail.com'
        receiver = 'gordon.s.kiesling@vanderbilt.edu'

        body = "There is a kickoff time for Western KY\nThe game will kick off at " + str(kickoff)

        # Structuring the email
        msg = MIMEText(body)
        msg['From'] = sender
        msg['To'] = receiver
        msg['Subject'] = "Kickoff time Announced !!!!!!!"

        # Access sender email and send email
        mail = smtplib.SMTP('smtp.gmail.com', 587)
        mail.ehlo()
        mail.starttls()
        mail.login('burner.checker@gmail.com', 'Rhsk24fn')
        message = msg.as_string()
        mail.sendmail(sender, receiver, message)
        mail.close()

    time.sleep(60 * 10)
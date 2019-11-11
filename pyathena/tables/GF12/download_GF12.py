import bs4
import requests
import os

if not os.path.isdir('./data'): os.mkdir('./data')

url = 'http://wise-obs.tau.ac.il/~orlyg/ion_by_ion/'
r = requests.get(url)
data = bs4.BeautifulSoup(r.text, "html.parser")
for l in data.find_all("a"):
    if l['href'].endswith('txt'):
        r = requests.get(url + l["href"])
        if r.status_code == 200:
            filename='./data/{}'.format(l["href"])
            with open(filename, 'wb') as f:
                f.write(r.content)
                print('Downloading... {} to {}'.format(l['href'],filename))

url = ['http://iopscience.iop.org/0067-0049/199/1/20/suppdata/apjs420150t2_ascii.txt',
       'http://wise-obs.tau.ac.il/~orlyg/cooling/CIEion/tab2.txt',
       'http://wise-obs.tau.ac.il/~orlyg/cooling/TDion/tab3.txt',
      ]
for url_ in url:
    r = requests.get(url_)
    if r.status_code == 200:
        filename='./data/'+os.path.basename(url_)
        with open(filename, 'wb') as f:
            f.write(r.content)
            print('Downloading... {}'.format(filename))

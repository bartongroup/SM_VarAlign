import urllib2
from retry import retry


# http://stackoverflow.com/questions/9446387/how-to-retry-urllib2-request-when-fails
@retry(urllib2.URLError, tries=4, delay=3, backoff=2)
def urlopen_with_retry(url):
    return urllib2.urlopen(url)
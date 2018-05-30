#!/usr/bin/python
import sys
#activate_this = '/var/www/reaction/python2/bin/activate_this.py'
#with open(activate_this) as file_:
#    exec(file_.read(), dict(__file__=activate_this))
#execfile(activate_this, dict(__file__=activate_this))
sys.path.insert(0,'/var/www/reaction/')
from init import app as application

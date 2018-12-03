#!/usr/bin/python

import sys
import SimpleHTTPServer
import SocketServer

class myHandler(SimpleHTTPServer.SimpleHTTPRequestHandler):
   def do_GET(self):
       self.send_response(301)
       self.send_header('Location','http://192.168.3.200:5000')
       self.end_headers()

if len(sys.argv) > 1:
   PORT = int( sys.argv[1] )
else:
   PORT = 5000

print "serving at port", PORT
handler = SocketServer.TCPServer(("", PORT), myHandler)
handler.serve_forever()

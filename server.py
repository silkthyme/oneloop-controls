import socket

# create the IP socket (endpoint in a communication)
# AF_INET == ipv4 (domain)
# SOCK_STREAM == TCP (TCP means connection-oriented)
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# bind to address
s.bind((socket.gethostname(), 1234))

# make a queue of 5 in case we have multiple incoming connections
s.listen(5)

while True:
	# now our endpoint knows about the other endpoint
	clientsocket, address = s.accept()
	print(f"Connection from {address} has been established.")


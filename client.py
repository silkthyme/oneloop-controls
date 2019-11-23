import socket

# make the TCP socket 
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# connect to server
s.connect((socket.gethostname(), 1234))

# receive data in a buffer size of 1024 bytes
message = s.recv(1024)

# print out data
print(message.decode("utf-8"))

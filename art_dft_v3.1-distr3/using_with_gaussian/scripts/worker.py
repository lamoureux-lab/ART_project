#!/usr/bin/python

import os
import sys
import time


class FileLock:
	def __init__(self, filename):
		self.filename = filename
		self.fd = None
		self.pid = os.getpid()

	def acquire(self):
		try:
			self.fd = os.open(self.filename, os.O_CREAT|os.O_EXCL|os.O_RDWR)
			# Only needed to let readers know who's locked the file
			os.write(self.fd, "%d" % self.pid)
			return 1    # return ints so this can be used in older Pythons
		except OSError:
			self.fd = None
			return 0

	def release(self):
		if not self.fd:
			return 0
		try:
			os.close(self.fd)
			os.remove(self.filename)
			return 1
		except OSError:
			return 0

	def __del__(self):
		self.release()


class PipedSet:

	def __init__(self, name = 'myfifo.dat'):
		self.lock = FileLock('.lock_'+name)
		self.name = name

		self.lock.acquire()
		if not os.path.exists(self.name):
			print 'Created file'
			f = open(self.name,'w')
			f.close()
		self.lock.release()
		
	def write(self, element):
		self.lock.acquire()
		print 'Writing file'
		with open(self.name, 'a+') as fout:
			for i in element:
				fout.write(i)
				time.sleep(0.1)
		self.lock.release()

	
	def read(self):
		self.lock.acquire()
		with open(self.name, 'r') as fin:
			out = fin.read()
		self.lock.release()
		return out


if __name__=='__main__':
	ps = PipedSet()
	for i in range(0,100):
		localList = ps.read()
		print 'Process %d read:'%os.getpid(), localList
		ps.write('This string should be intact')

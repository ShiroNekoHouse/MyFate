#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, random, base64, time
import ECCkeyGen as ECC
from federatedml.secureprotol import gmpy_math
import gmpy2
import test_pollard as pollard
import numpy as np
from federatedml.secureprotol.fixedpoint import FixedPointNumber

# p =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF 
# a =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC 
# b =0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93 
# n =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123 
# Gx=0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7 
# Gy=0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0
# K =100000000000000000#temp
# G =(Gx,Gy)

a	= 0
b	= 7
#modulo	= 2**256-2**32-2**9-2**8-2**7-2**6-2**4-1
p	= 115792089237316195423570985008687907853269984665640564039457584007908834671663
n	= 115792089237316195423570985008687907852837564279074904382605163141518161494337
#modulo	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
#order	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx	= 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy	= 32670510020758816978083085130507043184471273380659243275938904335757337482424
#Gx	= 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
#Gy	= 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8
K =100000000000000000#temp
G =(Gx,Gy)
PRECISION=10000000#temp

def encode(value):
	#encode a float to int,save the new int value and exponent
	if value==0:
		return (1,0)
	i=value
	expanse=1
	ans=0
	if i>0:
		while(i*10<1000000000):
			i=i*10.0
			expanse=expanse*10
		ans=int(round(i))
	else:
		i=-i
		while(i*10<1000000000):
			i=i*10
			expanse=expanse*10
		ans=int(round(i))
		expanse=-expanse
    #print(ans)
	return (expanse, ans)

def decode(encoding):
	value=encoding[1]
	expanse=encoding[0]
	ans=value/expanse
	return ans

def fate_encode(value,precision=None):
	if precision == None:
		return FixedPointNumber.encode(value,None,None,PRECISION)
	else:
		return FixedPointNumber.encode(value,None,None,precision)

def fate_decode(encode):
	return encode.decode()


def bytes2int(text,l,r): #加密输入l=31,解密输入l=32
	b=text
	data = []
	for i in range(r):
		a=b[:l]
		c = 0
		for j in range(l):
			c+=a[j]<<(8*(l-j-1))
		data.append(c)
		b=b[l:]
	return data
	
def int2bytes(data,l):    #解密输出l=31,加密输出l=32
	text = []
	for i in data:
		A = i
		for j in range(l)[::-1]:
			text.append((A >> 8*j) % 0x100)
	text = bytes(text)
	return text
	
def MessageDiv(A):
	if (len(A))%31!=0:
		A=A+bytes(31-(len(A))%31)
	r=len(A)//31
	return bytes2int(A,31,r)

class ECCKeypair(object):
	def __init__(self):
		pass

	@staticmethod
	def generate_keypair():
		"""return a new :class:`ECCPublicKey` and :class:`ECCPrivateKey`.
		"""
		PK, SK = generate_Keypair()  
		public_key=ECCPublicKey(PK)
		private_key=ECCPrivateKey(PK,SK)
		return public_key, private_key

class ECCPublicKey(object):
	"""Contains a public key and associated encryption methods.
	"""
	def __init__(self, publicKey):
		self.publicKey=publicKey

	def encrypt(self, value):
		"""Encode and Paillier encrypt a real number value.
		"""
		# encoding=encode(value)
		# exponent=encoding[0]
		# plaintext=encoding[1]
		
		#fate_encoding
		encoding=fate_encode(value)
		plaintext=encoding.encoding
		exponent=encoding.exponent
		print("plaintext:")
		print(plaintext)

		m_point=ECC.MultipyPoint(plaintext,G,a,p)
		flag = False
		while not flag:
			k=random.randrange(300,n-1)
			X1=ECC.MultipyPoint(k,G,a,p)
			X2=ECC.MultipyPoint(k,self.publicKey,a,p)
			if X2[0]!=None:
				flag=True
		
		# C=X2[0]*message%n
		C=ECC.PointAdd(a,p,m_point,X2)

		# data.append(X1[0])
		# data.append(X1[1])
		# data.append(C)
		data=ecc_number(X1,C,exponent,self)
		return data
		# return int2bytes(data,32)

class ECCPrivateKey(object):
	"""Contains a private key and associated decryption method.
	"""
	def __init__(self, publicKey, privateKey):
		self.publicKey=publicKey
		self.privateKey=privateKey
		
	def decrypt(self, encrypted_number):
		"""return the decrypted & decoded plaintext of encrypted_number.
		"""
		if encrypted_number==0:
			return 0
		# r=len(message)//96
		# C = bytes2int(message,32,r*3)
		X1=encrypted_number.X1
		C=encrypted_number.C
		# data = []
		
		# X1=(C[0],C[1])
		X2=ECC.MultipyPoint(self.privateKey,X1,a,p)
		C_y=(-C[1]+p)%p
		newC=(C[0],C_y)
		# V=ECC.modinv(X2[0], n)
		# # data.append((C[2]*V)%n)
		# plaintext=(C*V)%n
		tmpPoint=ECC.PointAdd(a,p,X2,newC)
		plaintext=pollard.pollard_method(tmpPoint)
		
		# plaintext=ECurvetoM(tmpPoint)
		#plaintext=decode((encrypted_number.exponent,plaintext))

		#fate_decode
		encoded=FixedPointNumber(plaintext,encrypted_number.exponent)
		plaintext=fate_decode(encoded)

		return plaintext
		# C=C[3:]
		# return int2bytes(data,31)


# def encrypt(m_point,PK):
# 	data = []
# 	flag = False
# 	while not flag:
# 		k=random.randrange(300,n-1)
# 		X1=ECC.MultipyPoint(k,G,a,p)
# 		X2=ECC.MultipyPoint(k,PK,a,p)
# 		if X2[0]!=None:
# 			flag=True
	
# 	# C=X2[0]*message%n
# 	C=ECC.PointAdd(a,p,m_point,X2)

# 	# data.append(X1[0])
# 	# data.append(X1[1])
# 	# data.append(C)
# 	data=ecc_number(X1,C)
# 	return data
# 	# return int2bytes(data,32)

# def decrypt(message,SK):
# 	# r=len(message)//96
# 	# C = bytes2int(message,32,r*3)
# 	X1=message.X1
# 	C=message.C
# 	# data = []
	
# 	# X1=(C[0],C[1])
# 	X2=ECC.MultipyPoint(SK,X1,a,p)
# 	C_y=(-C[1]+p)%p
# 	newC=(C[0],C_y)
# 	# V=ECC.modinv(X2[0], n)
# 	# # data.append((C[2]*V)%n)
# 	# plaintext=(C*V)%n
# 	tmpPoint=ECC.PointAdd(a,p,X2,newC)
# 	plaintext=pollard.pollard_method(tmpPoint)

# 	# plaintext=ECurvetoM(tmpPoint)
# 	return plaintext
# 	# C=C[3:]
# 	# return int2bytes(data,31)

def MtoECurve(m):
	j=random.randrange(0,K-1)
	x=m*K#+j
	if((m+1)*K>p):
		print("=================K is to big: error=================")
	A=(pow(x,3)+a*x+b)%p
	A_mpz=gmpy2.mpz(A)
	p_mpz=gmpy2.mpz(p)
	# print(gmpy_math.powmod(2,4,15))
	flag=gmpy2.powmod(A_mpz,(p_mpz-1)//2,p_mpz)   #pow(A,(n-1)/2)%n

	if(flag==1):#satisfy Quadratic residue 
		y=gmpy_math.powmod(A_mpz,(p_mpz+1)//4,p_mpz) #pow(A,(n+1)/4)%n
		return (x,y)
	else:
		print("=================not satisfy Quadratic residue================")

def ECurvetoM(m_point):
	return int(m_point[0]/K)

def sub_point(m_point):
	point_y=(-m_point[1]+p)%p
	new_point=(m_point[0],point_y)
	return new_point

def generate_Keypair():
	privateKey,publicKey=ECC.ECCKeyGen()
	return publicKey, privateKey

class ecc_number():
	def __init__(self,X1,C,exponent,publicKey):
		self.X1=X1
		self.C=C
		self.exponent=exponent
		self.publicKey=publicKey

	def increase_exponent_to(self, new_exponent):
		"""return ECCEncryptedNumber: 
			new ECCEncryptedNumber with same value but having great exponent.
		"""
		if new_exponent < self.exponent:
			raise ValueError("New exponent %i should be great than old exponent %i" % (new_exponent, self.exponent))
            
		factor = pow(FixedPointNumber.BASE, new_exponent - self.exponent)
		new_encryptednumber = self.__mul__(factor)        
		new_encryptednumber.exponent = new_exponent
		
		return new_encryptednumber
    
	def __align_exponent(self, x, y):
		"""return x,y with same exponet
		"""
		if x.exponent < y.exponent:
			x = x.increase_exponent_to(y.exponent)
		elif x.exponent > y.exponent:
			y = y.increase_exponent_to(x.exponent)
		
		return x, y
	
		
	# def increase_exponent_to(self, new_exponent):
	# 	"""return ecc_number: 
	# 		new ecc_number with same value but having great exponent.
	# 	"""

	# 	factor = int(abs(new_exponent)/abs(self.exponent))
	# 	new_encryptednumber = self.__mul__(factor) 
	# 	if(self.exponent<0):
	# 		new_encryptednumber.exponent = -abs(new_exponent)
	# 	else:
	# 		new_encryptednumber.exponent = abs(new_exponent)
			
	# 	return new_encryptednumber
	
	# def __align_exponent(self, x, y):
	# 	"""return x,y with same exponet
	# 	"""
	# 	if abs(x.exponent) < abs(y.exponent):
	# 		x = x.increase_exponent_to(y.exponent)
	# 	elif abs(x.exponent) > abs(y.exponent):
	# 		y = y.increase_exponent_to(x.exponent)
		
	# 	return x, y

	def __add__(self,other):
		if isinstance(other, ecc_number):
			return self.__add_encryptednumber(other)
		else:
			return self.__add_scalar(other)

	def __add_scalar(self,other):
		other=self.publicKey.encrypt(other)
		return self.__add_encryptednumber(other)

	def __add_encryptednumber(self,other):
		x, y = self.__align_exponent(self, other)
		# X1=ECC.PointAdd(a,p,self.X1,other.X1)
		# C=ECC.PointAdd(a,p,self.C,other.C)
		encryptednumber=self.__raw_add(x,y,x.exponent)
		return encryptednumber


	def __radd__(self, other):
		return self.__add__(other)

	def __raw_add(self,e_x,e_y,exponent):
		X1=ECC.PointAdd(a,p,e_x.X1,e_y.X1)
		C=ECC.PointAdd(a,p,e_x.C,e_y.C)
		return ecc_number(X1,C,exponent,self.publicKey)

	def __sub__(self, other):
		x, y = self.__align_exponent(self, other)
		# X1=ECC.PointAdd(a,p,self.X1,other.X1)
		# C=ECC.PointAdd(a,p,self.C,other.C)
		encryptednumber=self.__raw_sub(x,y,x.exponent)
		return encryptednumber
	
	def __raw_sub(self,e_x,e_y,exponent):
		print("????")
		C_y=(-e_y.C[1]+p)%p
		X1_y=(-e_y.X1[1]+p)%p
		newC=(e_y.C[0],C_y)
		newX1=(e_y.X1[0],X1_y)
		subEcc_num=ecc_number(newX1,newC,exponent,self.publicKey)
		return self.__add__(subEcc_num)

	def __rsub__(self, other):
		return other + (self * -1)

	def __rmul__(self,scalar):
		return self.__mul__(scalar)

	def __mul__(self,scalar):
		#change to int 
		if scalar == 0:
			return 0
		# if not isinstance(scalar, int):
		#point mul
		# encoding = encode(scalar)
		# expanse=encoding[0]
		# scalar=encoding[1]
		#fate_encode
		if scalar < 0:
			self = ecc_number(sub_point(self.X1),sub_point(self.C),self.exponent,self.publicKey)
			scalar = -scalar
		encoding=fate_encode(scalar,10)
		exponent=self.exponent+encoding.exponent
		scalar=encoding.encoding
		return ecc_number(ECC.MultipyPoint(scalar,self.X1,a,p),ECC.MultipyPoint(scalar,self.C,a,p),exponent,self.publicKey)

		# else:
		# 	if scalar>0:
		# 		flag=1
		# 	else:
		# 		flag=-1
		# 	return ecc_number(ECC.MultipyPoint(scalar,self.X1,a,p),ECC.MultipyPoint(scalar,self.C,a,p),self.exponent*flag,self.publicKey)




if __name__ == "__main__":
	publicKey,privateKey=ECCKeypair.generate_keypair()#SK,PK

	# choose=int(input('\n请选择明文编码:1.utf-8 2.GBK: '))-1
	# EnCo=coding[choose]
	# Message = bytes(input('请输入明文: '), encoding=EnCo)
	# m1=123456789
	# m2=123456789
	# mPoint=MtoECurve(m1)
	# mPoint2=MtoECurve(m2)

	# MALL=ECC.PointAdd(a,p,mPoint,mPoint2)
	# print(ECurvetoM(MALL))# bu manzu tongtai

	print("******************test****************")
	x_li = np.ones(2) * np.random.rand()
	print(x_li)
	aaa=0
	bbb=x_li[1]
	ccc=bbb-aaa

	# encoA=encode(aaa)
	# encoB=encode(bbb)
	# print (encoA)
	# print (encoB)
	# print(decode(encoA))
	# print(decode(encoB))

	# testx=12345678910
	# testy=12345678900
	# testPointX=ECC.MultipyPoint(aaa,G,a,p)
	# testPointY=ECC.MultipyPoint(bbb,G,a,p)
	# m_testPointX=(testPointX[0],(-testPointX[1]+p)%p)
	# testttt=ECC.MultipyPoint(bbb-aaa,G,a,p)
	# print(testttt)
	# pollard.pollard_method(testPointX)

	# testPointZ=ECC.MultipyPoint(testx+testy,G,a,p)
	# testPointX_Y=ECC.PointAdd(a,p,testPointY,m_testPointX)
	# print(testPointX_Y)
	print("**************************")

	# start = time.clock()
	ctext1 = publicKey.encrypt(aaa)
	ctext2 = publicKey.encrypt(bbb)
	print(ECC.PointAdd(a,p,ctext1.C,ctext2.C))

	# tmptext = base64.b64encode(ctext)
	# tmptext = str(tmptext, encoding="ascii")
	# print('\n密文(int形式输出):\n', ctext1)
	CTEXTALL = ctext1 - ctext2
	print(CTEXTALL.C)

	ptext = privateKey.decrypt(CTEXTALL)
	print('\n明文:\n')
	print(ptext)


	x_li = np.ones(10,dtype=int) * np.random.randint(100)
	y_li = np.ones(10,dtype=int) * np.random.randint(1000)        
	
	# for i in range(x_li.shape[0]):
	# 	x = x_li[i]
	# 	y = y_li[i]
	# 	testPointX=ECC.MultipyPoint(x,G,a,p)
	# 	testPointY=ECC.MultipyPoint(y,G,a,p)

	# 	en_x = encrypt(testPointX,publicKey)
	# 	en_y = encrypt(testPointY,publicKey)

		
	# 	en_res = en_x + en_y
		
	# 	res = x + y
		
	# 	de_en_res = decrypt(en_res,privateKey)
	# 	print(en_res)
	# 	print(res)

	# try:
	# 	plaintext = ptext.decode(coding[choose])
	# 	print(plaintext)
	# 	flag+=1
	# except:
	# 	#print("*"+ans[i]+"解码失败!")
	# 	pass
	# if flag==0:
	# 	print("\n解密失败!\n请核对密文/密钥的完整性或使用其他编码字符集\n")
	# elif flag==1:
	# 	print("")
	# else:
	# 	print("-----------------------------\n*请根据语义判断明文内容\n")


	# end = time.clock()
	# print("\n运算耗时 %f秒" % (end  - start))
	


###############################################################################
	# Inv = input('请选择:1.加密 2.解密: ')
	
	# if Inv == '2':
	# 	# d = eval(input('请输入私钥d: '))
	# 	d=k
	# 	Message = bytes(input('\n请输入base64格式的密文: '), encoding='ascii')
	# 	# start = time.clock()
	# 	Message = base64.b64decode(Message)
	# 	text = decrypt(Message,d)
	# else:
	# 	# Qx = eval(input('请输入公钥Qx: '))
	# 	# Qy = eval(input('请输入公钥Qy: '))
	# 	# Q = (Qx,Qy)
	# 	Q=D
	# 	choose=int(input('\n请选择明文编码:1.utf-8 2.GBK: '))-1
	# 	EnCo=coding[choose]
	# 	Message = bytes(input('请输入明文: '), encoding=EnCo)
	# 	# start = time.clock()
	# 	text = encrypt(Message,Q[0],Q[1])
	
	# ans = ['UTF-8','GBK']
	
	# if Inv == '1':
	# 	text = base64.b64encode(text)
	# 	text = str(text, encoding="ascii")
	# 	print('\n密文(以base64形式输出):\n', text)
	# else:
	# 	#print(text)
	# 	#print(len(text))
	# 	flag=0
	# 	print('\n明文:\n')
	# 	for i in range(2):
	# 		try:
	# 			plaintext = " "+ans[i]+': '+text.decode(coding[i])
	# 			print(plaintext)
	# 			flag+=1
	# 		except:
	# 			#print("*"+ans[i]+"解码失败!")
	# 			pass
	# 	if flag==0:
	# 		print("\n解密失败!\n请核对密文/密钥的完整性或使用其他编码字符集\n")
	# 	elif flag==1:
	# 		print("")
	# 	else:
	# 		print("-----------------------------\n*请根据语义判断明文内容\n")
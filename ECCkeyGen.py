#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os, random, base64, time

def egcd(a, b):
	if a == 0:
		return (b, 0, 1)
	else:
		g, y, x = egcd(b % a, a)  # g为公因子
		return (g, x - (b // a) * y, y)

def modinv(a, m):
	g, x, y = egcd(a, m)
	if (g != 1):
		print('modular inverse does not exist')
	else:
		return x % m

def PointAdd(a,p,A,B):
	if A[0]==-1:
		return B
	if B[0]==-1:
		return A
	if A[0]==B[0]:
		if A[1]!=B[1]:
			return(-1,-1)
		else:
			lam=(((3*(A[0]**2)+a)%p)* modinv(2*A[1], p))%p
	else:
		lam=(((B[1]-A[1])%p)* modinv((B[0]-A[0])%p, p))%p
	x3=(lam**2-A[0]-B[0] )% p
	y3=(lam*(A[0]-x3)-A[1]) % p
	return (x3,y3)
	
def MultipyPoint(n,A,a,p):   #借鉴了模重复平方计算法
	D=(-1,-1)
	E=bin(n)[2:]
	for i in range(len(E)):
		D=PointAdd(a,p,D,D)
		if E[i]=="1":
			D=PointAdd(a,p,D,A)
	return D
	
# p =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF 
# a =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC 
# b =0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93 
# n =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123 
# Gx=0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7 
# Gy=0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0
# K =100000000000000000#temp
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
# G =(Gx,Gy)

def ECCKeyGen():
	k=random.randrange(300,n-1)
	G=(Gx,Gy)
	D=MultipyPoint(k,G,a,p)
	return k,D

if __name__ == "__main__":

	start = time.clock()
	
	k=random.randrange(300,n-1)
	A=(Gx,Gy)
	D=MultipyPoint(k,A,a,p)
	print(D)
	#print(hex(k).upper(),D)
	
	end = time.clock()
	print("公钥:Qx =","0x"+hex(D[0]).upper()[2:],"\n     Qy =","0x"+hex(D[1]).upper()[2:])
	print("私钥: d =","0x"+hex(k).upper()[2:])
	print("\n运算耗时 %f秒" % (end  - start))
	os.system("PAUSE")
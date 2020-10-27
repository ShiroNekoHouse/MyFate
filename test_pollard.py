import pollard_kangaroo as pollard
import gmpy2

pow2bits	= 32	# bits/suborder/exp key

pubkeys = {
	  16: ('029d8c5d35231d75eb87fd2c5f05f65281ed9573dc41853288c62ee94eb2590b7a', 0xc936)
	, 24: ('036ea839d22847ee1dce3bfc5b11f6cf785b0682db58c35b63d1342eb221c3490c', 0xdc2a04)
	, 32: ('0209c58240e50e3ba3f833c82655e8725c037a2294e14cf5d73a5df8d56159de69', 0xb862a62e)
	, 33: ('02ed949eaca31df5e8be9bf46adc1dfae1734b8900dcc303606831372955c728da', False) #0x01abcd1234
	, 40: ('03a2efa402fd5268400c77c20e574ba86409ededee7c4020e4b9f0edbee53de0d4', 0xe9ae4933d6)
	, 50: ('03f46f41027bbf44fafd6b059091b900dad41e6845b2241dc3254c7cdd3c5a16c6', 0x022bd43c2e9354)
	, 55: ('0385a30d8413af4f8f9e6312400f2d194fe14f02e719b24c3f83bf1fd233a8f963', 0x6abe1f9b67e114)
	, 60: ('0348e843dc5b1bd246e6309b4924b81543d02b16c8083df973a89ce2c7eb89a10d', 0x0FC07A1825367BBE)
	, 70: ('0290e6900a58d33393bc1097b5aed31f2e4e7cbd3e5466af958665bc0121248483', 0x349B84B6431A6C4EF1)
	, 80: ('037e1238f7b1ce757df94faa9a2eb261bf0aeb9f84dbf81212104e78931c2a19dc', 0xEA1A5C66DCC11B5AD180)
	, 90: ('035c38bd9ae4b10e8a250857006f3cfd98ab15a6196d9f4dfd25bc7ecc77d788d5', 0x02CE00BB2136A445C71E85BF)
	,100: ('03d2063d40402f030d4cc71331468827aa41a8a09bd6fd801ba77fb64f8e67e617', 0x0af55fc59c335c8ec67ed24826)
	,105: ('03bcf7ce887ffca5e62c9cabbdb7ffa71dc183c52c04ff4ee5ee82e0c55c39d77b', False)
}

# p =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF 
# a =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC 
# b =0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93 
# n =0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123 
# Gx=0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7 
# Gy=0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0

# A_curve	= gmpy2.mpz(a)
# B_curve	= gmpy2.mpz(b)
# modulo	= gmpy2.mpz(p)
# order	= gmpy2.mpz(n)
# Gx	= gmpy2.mpz(Gx)
# Gy	= gmpy2.mpz(Gy)

# A_curve	= 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC 
# B_curve	= 0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93 
# #modulo	= 2**256-2**32-2**9-2**8-2**7-2**6-2**4-1
# modulo	= 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
# order	= 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123
# #modulo	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
# #order	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
# Gx	= 0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7
# Gy	= 0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0
# #Gx	= 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
# #Gy	= 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8

A_curve	= 0
B_curve	= 7
#modulo	= 2**256-2**32-2**9-2**8-2**7-2**6-2**4-1
modulo	= 115792089237316195423570985008687907853269984665640564039457584007908834671663
order	= 115792089237316195423570985008687907852837564279074904382605163141518161494337
#modulo	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
#order	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx	= 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy	= 32670510020758816978083085130507043184471273380659243275938904335757337482424
#Gx	= 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
#Gy	= 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8


A_curve	= gmpy2.mpz(A_curve)
B_curve	= gmpy2.mpz(B_curve)
modulo	= gmpy2.mpz(modulo)
order	= gmpy2.mpz(order)
Gx	= gmpy2.mpz(Gx)
Gy	= gmpy2.mpz(Gy)

def pollard_method(pubkey0):

	# python2+gmpy2 speed-up +8%


	Gp = pollard.Point(Gx,Gy) 
	Zp = pollard.Point(0,0)
	Sp = [Gp]
	for k in pollard.xrange(255): 
		Sp.append(pollard.mul_2a(Sp[k]))

	flag_pow2bits = True
	flag_keyspace = False

	# pubkey0, prvkey0 = pubkeys[pow2bits]
	# pubkey0 = pollard.mul_ka(prvkey0,Gp)

	# check format pubkey
	# if len(pubkey0)==130:
	#     X = int(pubkey0[2:66], 16)
	#     Y = int(pubkey0[66:],16)
	#     flag_compress = False
	#     print("[format] uncompressed")
	# elif len(pubkey0)==66:
	#     X = int(pubkey0[2:66], 16)
	#     # calculation Y from X if pubkey is compressed
	#     Y = pollard.getX2Y(X,int(pubkey0[:2])-2)
	#     flag_compress = True
	#     print("[format] compressed")
	# else:
	#     print("[error] pubkey len(66/130) invalid!")
	#     # usage()
	X=pubkey0[0]
	Y=pubkey0[1]

	#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
	print("[Xcoordinate] %064x" % X)
	print("[Ycoordinate] %064x" % Y)
	#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

	if not pollard.is_on_curve(X,Y):
		print("[error] the given point not lies on the elliptic curve!")
		pollard.usage()
	else:
		print("they are on the curve")
	# wild root
	W0p = pollard.Point(X,Y)

	L = 2**(pow2bits-1)
	U = 2**pow2bits

	W = U - L
	Wsqrt = W**0.5
	#Wsqrt = math.sqrt(W)
	Wsqrt = int(Wsqrt)


	# M == (L+U)/2 == L+(W/2)
	#M = (L + U)//2
	M = L + (W//2)

	pow2L = pow2bits-1
	pow2U = pow2bits
	pow2W = pow2bits-1

	prvkey, runjump, runtime, lenT,lenW = pollard.KANGAROOS(Gp,Sp,Wsqrt,pow2W,M,W0p)

	return prvkey
"""Paillier encryption library for partially homomorphic encryption."""

#
#  Copyright 2019 The FATE Authors. All Rights Reserved.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

from collections.abc import Mapping
from federatedml.secureprotol.fixedpoint import FixedPointNumber
from federatedml.secureprotol import gmpy_math
import random
import binascii
from random import choice
import sm3, ecc_func
# 选择素域，设置椭圆曲线参数

default_ecc_table = {
    'n': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123',
    'p': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF',
    'g': '32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7'\
         'bc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0',
    'a': 'FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC',
    'b': '28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93',
}
PRIVATE_KEY = '00B9AB0B828FF68872F21A837FC303668428DEA11DCD1B24429D0C99E24EED83D5'
PUBLIC_KEY = 'B9C9A6E04E9C91F7BA880429273747D7EF5DDEB0BB2FF6317EB00BEF331A83081A6994B8993F3F5D6EADDDB81872266C87C018FB4162F5AF347B483E24620207'

class CryptSM2(object):
    def __init__(self, private_key=None, public_key=None, ecc_table=default_ecc_table):
        private_key = '00B9AB0B828FF68872F21A837FC303668428DEA11DCD1B24429D0C99E24EED83D5'
        public_key = 'B9C9A6E04E9C91F7BA880429273747D7EF5DDEB0BB2FF6317EB00BEF331A83081A6994B8993F3F5D6EADDDB81872266C87C018FB4162F5AF347B483E24620207'
        self.private_key = private_key
        self.public_key = public_key
        self.para_len = len(ecc_table['n'])
        self.ecc_a3 = (
            int(ecc_table['a'], base=16) + 3) % int(ecc_table['p'], base=16)
        self.ecc_table = ecc_table

    
    @staticmethod
    def generate_keypair(n_length=1024):
        """return a new :class:`ECCPublicKey` and :class:`ECCPrivateKey`.
        call EC_KEY_new()
        """
        pass 

    def _kg(self, k, Point):  # kP运算
        Point = '%s%s' % (Point, '1')
        mask_str = '8'
        for i in range(self.para_len - 1):
            mask_str += '0'
        mask = int(mask_str, 16)
        Temp = Point
        flag = False
        for n in range(self.para_len * 4):
            if (flag):
                Temp = self._double_point(Temp)
            if (k & mask) != 0:
                if (flag):
                    Temp = self._add_point(Temp, Point)
                else:
                    flag = True
                    Temp = Point
            k = k << 1
        return self._convert_jacb_to_nor(Temp)

    def _double_point(self, Point):  # 倍点
        l = len(Point)
        len_2 = 2 * self.para_len
        if l< self.para_len * 2:
            return None
        else:
            x1 = int(Point[0:self.para_len], 16)
            y1 = int(Point[self.para_len:len_2], 16)
            if l == len_2:
                z1 = 1
            else:
                z1 = int(Point[len_2:], 16)

            T6 = (z1 * z1) % int(self.ecc_table['p'], base=16)
            T2 = (y1 * y1) % int(self.ecc_table['p'], base=16)
            T3 = (x1 + T6) % int(self.ecc_table['p'], base=16)
            T4 = (x1 - T6) % int(self.ecc_table['p'], base=16)
            T1 = (T3 * T4) % int(self.ecc_table['p'], base=16)
            T3 = (y1 * z1) % int(self.ecc_table['p'], base=16)
            T4 = (T2 * 8) % int(self.ecc_table['p'], base=16)
            T5 = (x1 * T4) % int(self.ecc_table['p'], base=16)
            T1 = (T1 * 3) % int(self.ecc_table['p'], base=16)
            T6 = (T6 * T6) % int(self.ecc_table['p'], base=16)
            T6 = (self.ecc_a3 * T6) % int(self.ecc_table['p'], base=16)
            T1 = (T1 + T6) % int(self.ecc_table['p'], base=16)
            z3 = (T3 + T3) % int(self.ecc_table['p'], base=16)
            T3 = (T1 * T1) % int(self.ecc_table['p'], base=16)
            T2 = (T2 * T4) % int(self.ecc_table['p'], base=16)
            x3 = (T3 - T5) % int(self.ecc_table['p'], base=16)

            if (T5 % 2) == 1:
                T4 = (T5 + ((T5 + int(self.ecc_table['p'], base=16)) >> 1) - T3) % int(self.ecc_table['p'], base=16)
            else:
                T4 = (T5 + (T5 >> 1) - T3) % int(self.ecc_table['p'], base=16)

            T1 = (T1 * T4) % int(self.ecc_table['p'], base=16)
            y3 = (T1 - T2) % int(self.ecc_table['p'], base=16)

            form = '%%0%dx' % self.para_len
            form = form * 3
            return form % (x3, y3, z3)

    def _add_point(self, P1, P2):  # 点加函数，P2点为仿射坐标即z=1，P1为Jacobian加重射影坐标
        len_2 = 2 * self.para_len
        l1 = len(P1)
        l2 = len(P2)
        if (l1 < len_2) or (l2 < len_2):
            return None
        else:
            X1 = int(P1[0:self.para_len], 16)
            Y1 = int(P1[self.para_len:len_2], 16)
            if (l1 == len_2):
                Z1 = 1
            else:
                Z1 = int(P1[len_2:], 16)
            x2 = int(P2[0:self.para_len], 16)
            y2 = int(P2[self.para_len:len_2], 16)

            T1 = (Z1 * Z1) % int(self.ecc_table['p'], base=16)
            T2 = (y2 * Z1) % int(self.ecc_table['p'], base=16)
            T3 = (x2 * T1) % int(self.ecc_table['p'], base=16)
            T1 = (T1 * T2) % int(self.ecc_table['p'], base=16)
            T2 = (T3 - X1) % int(self.ecc_table['p'], base=16)
            T3 = (T3 + X1) % int(self.ecc_table['p'], base=16)
            T4 = (T2 * T2) % int(self.ecc_table['p'], base=16)
            T1 = (T1 - Y1) % int(self.ecc_table['p'], base=16)
            Z3 = (Z1 * T2) % int(self.ecc_table['p'], base=16)
            T2 = (T2 * T4) % int(self.ecc_table['p'], base=16)
            T3 = (T3 * T4) % int(self.ecc_table['p'], base=16)
            T5 = (T1 * T1) % int(self.ecc_table['p'], base=16)
            T4 = (X1 * T4) % int(self.ecc_table['p'], base=16)
            X3 = (T5 - T3) % int(self.ecc_table['p'], base=16)
            T2 = (Y1 * T2) % int(self.ecc_table['p'], base=16)
            T3 = (T4 - X3) % int(self.ecc_table['p'], base=16)
            T1 = (T1 * T3) % int(self.ecc_table['p'], base=16)
            Y3 = (T1 - T2) % int(self.ecc_table['p'], base=16)

            form = '%%0%dx' % self.para_len
            form = form * 3
            return form % (X3, Y3, Z3)

    def _convert_jacb_to_nor(self, Point): # Jacobian加重射影坐标转换成仿射坐标
        len_2 = 2 * self.para_len
        x = int(Point[0:self.para_len], 16)
        y = int(Point[self.para_len:len_2], 16)
        z = int(Point[len_2:], 16)
        z_inv = pow(z, int(self.ecc_table['p'], base=16) - 2, int(self.ecc_table['p'], base=16))
        z_invSquar = (z_inv * z_inv) % int(self.ecc_table['p'], base=16)
        z_invQube = (z_invSquar * z_inv) % int(self.ecc_table['p'], base=16)
        x_new = (x * z_invSquar) % int(self.ecc_table['p'], base=16)
        y_new = (y * z_invQube) % int(self.ecc_table['p'], base=16)
        z_new = (z * z_inv) % int(self.ecc_table['p'], base=16)
        if z_new == 1:
            form = '%%0%dx' % self.para_len
            form = form * 2
            return form % (x_new, y_new)
        else:
            return None

    def encrypt(self, m):
        # 加密函数，data消息(bytes)
        # msg = data.hex() # 消息转化为16进制字符串
        k = func.random_hex(self.para_len)
        M=self._kg(m,self.ecc_table['g'])
        C1 = self._kg(int(k,16),self.ecc_table['g'])
        kP = self._kg(int(k,16),self.public_key)
        C2=self._add_point(M,kP)
        return [C1,C2]
        # return bytes.fromhex('%s%s' % (C1,C2))


    def decrypt(self, c):
        # 解密函数，data密文（bytes）
        # data = data.hex()
        if not isinstance(encrypted_number, ECCEncryptedNumber):
            raise TypeError("encrypted_number should be an PaillierEncryptedNumber, \
                             not: %s" % type(encrypted_number))

        # len_2 = 2 * self.para_len
        # len_3 = len_2 + 64
        # C1 = c[0:len_2]
        # # C3 = data[len_2:len_3]
        # C2 = c[len_2:]
        C1=c.C1
        C2=c.C2
        xy = self._kg(int(self.private_key,16),C1)

        M=self._add_point(C2,-xy)
        return bytes.fromhex(M)


class ECCPublicKey(object):
    """Contains a public key and associated encryption methods.
    """
    def __init__(self, public_key):
        self.sm2_crypt=CryptSM2()
        self.public_key=PUBLIC_KEY

    
    def encrypt(self, value):
        """Encode and ECC encrypt a real number value.
        """
        ciphertext=self.sm2_crypt.encrypt(value)
        encryptednumber = ECCEncryptedNumber(self.public_key, ciphertext[0],ciphertext[1])
        return encryptednumber
   

class ECCPrivateKey(object):
    """Contains a private key and associated decryption method.
    """
    def __init__(self,private_key):
        self.sm2_crypt = CryptSM2()
        self.private_key=PRIVATE_KEY

   
    def decrypt(self, encrypted_number):
        """return the decrypted & decoded plaintext of encrypted_number.
        """        
        if not isinstance(encrypted_number, ECCEncryptedNumber):
            raise TypeError("encrypted_number should be an PaillierEncryptedNumber, \
                             not: %s" % type(encrypted_number))

        decrypt_value=self.sm2_crypt.decrypt(encrypted_number)
        
        return decrypt_value

class ECCEncryptedNumber(object):
    """Represents the Paillier encryption of a float or int.
    """
    def __init__(self, public_key, C1, C2):
        self.public_key = public_key
        self.C1 = C1
        self.C2 = C2
        self.sm2_crypt = CryptSM2()

    def ciphertext(self, be_secure=True):
        """return the ciphertext of the ECCEncryptedNumber.
        """

        return [self.C1, self.C2]

  
    def __add__(self, other):       
        if isinstance(other, ECCEncryptedNumber):
            return self.__add_encryptednumber(other)
        else:
            return self.__add_scalar(other)

    def __radd__(self, other):        
        return self.__add__(other)

    def __sub__(self, other):
        return self + (other * -1)

    def __rsub__(self, other):
        return other + (self * -1)
    
    def __rmul__(self, scalar):
        return self.__mul__(scalar)
    
    def __truediv__(self, scalar):
        return self.__mul__(1 / scalar)  
    
    def __mul__(self, scalar):
        """return Multiply by an scalar(such as int, float)
        """
        C1=sm2_crypt._kg(int(scalar,self.C1)
        C2=sm2_crypt._kg(int(scalar,self.C2)

        return ECCEncryptedNumber(self.public_key, C1,C2)  
    
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
    
    def __add_scalar(self, scalar):
        """return ECCEncryptedNumber: z = E(x) + y
        """
        encoded = FixedPointNumber.encode(scalar, 
                                          self.public_key.n, 
                                          self.public_key.max_int,                                          
                                          max_exponent=self.exponent)

        return self.__add_fixpointnumber(encoded)     
        
    def __add_fixpointnumber(self, encoded):
        """return ECCEncryptedNumber: z = E(x) + FixedPointNumber(y)
        """
        if self.public_key.n != encoded.n:
            raise ValueError("Attempted to add numbers encoded against different public keys!")

        # their exponents must match, and align.
        x, y = self.__align_exponent(self, encoded)

        encrypted_scalar = x.public_key.raw_encrypt(y.encoding, 1)
        encryptednumber = self.__raw_add(x.ciphertext(False), encrypted_scalar, x.exponent)
        
        return encryptednumber

    def __add_encryptednumber(self, other):
        """return ECCEncryptedNumber: z = E(x) + E(y)
        """
        C_x_plus_y_1=self.sm2_crypt._add_point(e_x.C1,e_y.C1)
        C_x_plus_y_2=self.sm2_crypt._add_point(e_x.C2,e_y.C2)
        
        return ECCEncryptedNumber(self.public_key, ciphertext, exponent)
        
        return encryptednumber

    def __raw_add(self, e_x, e_y, exponent):
        """return the integer E(x + y) given ints E(x) and E(y).
        """
        C_x_plus_y_1=self.sm2_crypt._add_point(e_x.C1,e_y.C1)
        C_x_plus_y_2=self.sm2_crypt._add_point(e_x.C2,e_y.C2)
        
        return ECCEncryptedNumber(self.public_key, ciphertext, exponent)
        
        return encryptednumber

        
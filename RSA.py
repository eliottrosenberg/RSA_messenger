## Simple RSA Messenger by Eliott Rosenberg, 2021.
## See https://en.wikipedia.org/wiki/RSA_(cryptosystem)#Operation

## requirements:
# sympy

import math
import sys
from typing import Tuple
sys.setrecursionlimit(1000000)


class sender:
    def __init__(self,partner_n=None,partner_e=None):
        
        self.partner_n = partner_n
        self.partner_e = partner_e
    
    def encrypt(self,message,unicode=False):
        if self.partner_n == None or self.partner_e == None:
            raise Exception('Need public key from partner.')
            return False
        message_number = string_to_number(message,unicode)
        c_to_send = pow(message_number,self.partner_e,self.partner_n)
        return c_to_send
    

class receiver:
    def __init__(self,n_bits=2048,my_n=None,my_e=None,my_d=None,generate_key=True):
        if (my_n == None or my_e == None or my_d == None) and generate_key:
            self.set_params(n_bits)
        else:
            self.my_n = my_n
            self.my_e = my_e
            self.my_d = my_d
            
    def set_params(self,n_bits):
        p, q = pick_primes(n_bits)
        self.set_params_from_primes(p,q)
    
    def set_params_from_primes(self,p,q):
        self.my_n = p*q
        lam = (p-1)*(q-1)//math.gcd(p-1,q-1)
        self.my_e = 65537
        if self.my_e > lam:
            while not (self.my_e < lam and math.gcd(self.my_e,lam) == 1):
                self.my_e -= 1
        while not (self.my_e < lam and math.gcd(self.my_e,lam) == 1):
            self.my_e += 1
        self.my_d = modinv(self.my_e, lam)  # in python 3.8+: self.my_d = pow(self.my_e,-1,lam)
        
    def decrypt(self,c_received,unicode=False):
        message_number = pow(c_received,self.my_d,self.my_n)
        message = number_to_string(message_number,unicode)
        return message
        
    def generate_sender(self):
        return sender(self.my_n,self.my_e)


def pick_primes(n_bits):
    # picks primes with a product near the desired number of bits
    from sympy.ntheory.generate import randprime
    p = randprime(2**(n_bits//2),2**(n_bits//2 +1))
    q = randprime(2**(n_bits//2),2**(n_bits//2 +1))
    mult = 1.1
    while q == p:
        q = randprime(2**(n_bits//2),mult*2**(n_bits//2 +1))
        mult = mult*1.1
    return p, q
    
def string_to_number(string,unicode=False):
    if unicode:
        bits_per_block = 18
    else:
        bits_per_block = 7
    message_binary = ''
    for c in string:
        c_binary = bin(ord(c))[2:]
        n_zeros = bits_per_block - len(c_binary)
        c_binary = ''.join(['0' for _ in range(n_zeros)]) + c_binary
        message_binary += c_binary
    return int(message_binary,2)
    
def number_to_string(number,unicode=False):
    if unicode:
        bits_per_block = 18
    else:
        bits_per_block = 7
    message_string = ''
    number_binary = bin(number)[2:]
    extra_zeros = (bits_per_block - len(number_binary)%bits_per_block)%bits_per_block
    number_binary = ''.join(['0' for _ in range(extra_zeros)]) + number_binary
    num_characters = len(number_binary)//bits_per_block
    
    string = ''
    for which_character in range(num_characters):
        character_bin = number_binary[(which_character*bits_per_block):( (which_character+1)*bits_per_block)]
        string += chr(int(character_bin,2))
        
    return string
    

# the following are from https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Extended_Euclidean_algorithm

def egcd(a: int, b: int) -> Tuple[int, int, int]:
    ## return (g, x, y) such that a*x + b*y = g = gcd(a, b)
    if a == 0:
        return (b, 0, 1)
    else:
        b_div_a, b_mod_a = divmod(b, a)
        g, x, y = egcd(b_mod_a, a)
        return (g, y - b_div_a * x, x)

def modinv(a: int, b: int) -> int:
    ## return x such that (x * a) % b == 1
    g, x, _ = egcd(a, b)
    if g != 1:
        raise Exception('gcd(a, b) != 1')
    return x % b
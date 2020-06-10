#!/bin/env python


"""Extract a sub PSSM matching a fragment of a protein structure
"""

__author__ = 'Rosalie Lipsh-Sokolik'


import argparse
import os
import pandas as pd
import re
import logging
import os
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import PPBuilder


def get_pdb_seq(path):
    """
    :param path: path to a pdb file
    """
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            file_content = f.read().decode('utf-8')
            tmp = 'tmp'
            open(tmp, 'w').write(file_content)
            pdb = PDBParser().get_structure('', tmp)
            os.remove(tmp)
    else:
        pdb = PDBParser().get_structure('', path)
    pp = PPBuilder().build_peptides(pdb)
    seq = ''.join([str(i.get_sequence()) for i in pp])
    return seq


class PSSM:                                                                          
    cols = ['resseq', 'aa', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
            'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
            'A_cons', 'R_cons', 'N_cons', 'D_cons', 'C_cons', 'Q_cons', 'E_cons', 
            'G_cons', 'H_cons', 'I_cons', 'L_cons', 'K_cons', 'M_cons', 'F_cons', 
            'P_cons', 'S_cons', 'T_cons', 'W_cons', 'Y_cons', 'V_cons', 
            'num1', 'num2']

    def __init__(self, path=''):
        self.header = []
        self.footer = []
        self.pssm = None
        if path:
            self.parse_pssm(path)
        
    def parse_pssm(self, pssm_path):
        """Parses a PSSM file from """
        pssm = open(pssm_path, 'r').readlines()
        self.pssm = pd.DataFrame(columns=PSSM.cols)
        flag_header = True
        for line in pssm:
            if self.is_pssm_line(line):
                flag_header = False
                self.pssm.loc[len(self.pssm)] = self.__parse_pssm_line(line)
            elif flag_header:
                self.header.append(line)
            else:
                self.footer.append(line)
    
    def __parse_pssm_line(self, line):
        """
        If badly formatted, assumes the pseudocount are connected for some
        reason
        """
        result = line.split()
        if len(result) == 43:  # assuming messed up pseudocounts
            result = result[:-1] + [result[-1][:-5], result[-1][-5:]]
        if len(result) != 44:
            raise ValueError('Incorrect PSSM line:\n{}'.format(line))
        return result

    @staticmethod
    def is_pssm_line(line):  
        """A line is not a PSSM line if:
            1. has 2 letters in a row
            2. empty
            3. has no numbers
        """
        return not (re.search('[a-zA-Z]{2,}', line) or
                    line == '\n' or
                    not re.search('\d', line))

    def __get_lines(self, start=None, end=None):
        """
        Indexing as python lists: starts from 0 and [start:end).
        :return: dataframe with wanted lines
        """
        return self.pssm[start: end]

    def get_lines(self, start=None, end=None):
        """
        :param start and end:'starts from 0 and [start:end)
        :return: a new PSSM object
        """
        new_pssm = PSSM()
        new_pssm.header = self.header
        new_pssm.footer = self.footer
        new_pssm.pssm = self.__get_lines(start, end)
        return new_pssm

    def get_sequence(self):
        """Returns the aa column as a string"""
        return ''.join(self.pssm.aa.tolist())

    def write_pssm(self, path):
        pssm_str = '{:>5} {}  ' + '{:>4}'*20 + ' ' + '{:>4}'*20 + '{:>6}{:>5}\n'
        f = open(path, 'w')
        f.writelines(self.header)
        for index, data in self.pssm.iterrows():
            f.write(pssm_str.format(*data.tolist()))
        f.writelines(self.footer)

    def cut_pssm_by_seq(self, seq):
        """Cuts the pssms according to the sequence in the pdb
        :param seq: a sequence of a sub pdb of the original structure of the 
        pssm. if it is a file, assumes it is a pdb structure and reads its
        sequence.
        :return: the sub pssm corresponding for the pdb fragment
        """
        seq = get_pdb_seq(seq)
        pssm_seq = self.get_sequence()
        count = pssm_seq.count(seq)
        if count > 1:
            msg = 'The seq appears {} times in the pssm:\n{}\n{}'.format(
                count, seq, pssm_seq)
            raise Exception(msg)
        elif count < 1:
            msg = 'The seq is not in the pssm:\n{}\n{}'.format(seq, pssm_seq)
            raise Exception(msg)
        startI = pssm_seq.index(seq)
        endI = startI + len(seq)
        pssm = self.get_lines(startI, endI)
        assert pssm.get_sequence() == seq
        return pssm
    
    
def parse_args():                                                                                               
    """"""
    desc = 'Cuts the blades\'s lines from original pssm (of the full protein).'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('pdb', help='Path to truncated pdb. ')
    parser.add_argument('pssm', help='Path to pssm of full protein')
    parser.add_argument('-name', default='', help='name for the new pssm file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    pssm_org = PSSM(path=args.pssm)
    pssm = pssm_org.cut_pssm_by_seq(args.pdb)
    print(pssm.pssm)
    name = args.name if args.name else os.path.basename(args.pdb).replace('.gs', '').replace('.pdb', '.pssm')
    pssm.write_pssm(name)
    print(name)

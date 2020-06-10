"""
Refining fragments to align better to template
"""


#!/usr/bin/env python3

import os
import logging
import argparse
import sys
from flab.utils.parse_actions import FullPaths, FullPathsList


__author__ = 'Gideon Lapidoth'
__maintainer__ = 'Rosalie Lipsh'


LGR = logging.getLogger(__file__)


class BestMatchFinder:
    """
    """
    def __init__(self, template, pdbs, save_path=''):
        cmd.delete('all')
        self.template = template
        self.pdbs = pdbs if isinstance(pdbs, list) else [pdbs]
        
        self.templateN = 'template'
        self.targetN = 'target'
        self.template_start_strand = 'template_start_strand'
        self.template_end_strand = 'template_end_strand'
        self.target_start_strand = 'target_start_strand'
        self.target_end_strand = 'target_end_strand'

        self.__set_template()
        self.save_path = save_path
        if not save_path:
            self.save_path = os.path.join(os.getcwd(), 'pdbs')
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

    def __set_template(self):
        """Loading the template structure, and define the start and end of the
        template strand by first and last 3 residues
        """
        cmd.load(self.template, object=self.templateN)
        space = {'residues': list()}
        cmd.iterate(self.templateN, 'residues.append(resi)', space=space)
        residues = sorted([int(r) for r in set(space['residues'])])
        selection = 'resi {}-{} and {} and name CA+C+O+N'
        cmd.select(self.template_start_strand,
                   selection.format(residues[0], residues[2], self.templateN))
        cmd.select(self.template_end_strand,
                   selection.format(residues[-3], residues[-1],
                                    self.templateN))

    def __is_continues(self, i, residues):
        """Makes sure a strech of 3 resiues is continues"""
        return ((residues[i+1] - residues[i] == 1) and
                (residues[i+2] - residues[i] == 2))
    
    def __pair_fit(self, start1, start2, end1, end2):
        try:
            selection = '{} and resi {}-{} and name c+ca+n+o'
            cmd.select(self.target_start_strand, selection.format(
                self.targetN, start1, start2))
            cmd.select(self.target_end_strand, selection.format(
                self.targetN, end1, end2))
        except:
            return -1
        target_pair = '{} + {}'.format(self.target_start_strand,
                                       self.target_end_strand)
        template_pair = '{} + {}'.format(self.template_start_strand,
                                         self.template_end_strand)
        score = cmd.pair_fit(target_pair, template_pair)
        return score
    
    def __find_best_match(self, pdb, high_threshold=2, prefix='fbm'):
        best = 10
        bestI = {'start1': -1, 'start2': -1, 'end1': -1, 'end2': -1}
        best_start = -1
        best_end = -1
        cmd.load(pdb, object=self.targetN)
        cmd.remove("het")
        space = {'residues': list()}
        cmd.iterate('{} and not ss h'.format(self.targetN),
                    'residues.append(resi)', space=space)
        target_strand = sorted([int(r) for r in set(space['residues'])])
        for j in range(0, len(target_strand)-5):
            if not self.__is_continues(j, target_strand):
                continue
            for k in range(j+3, len(target_strand)-2):
                if not self.__is_continues(k, target_strand):
                    continue
                score = self.__pair_fit(target_strand[j], target_strand[j+2],
                                        target_strand[k], target_strand[k+2])
                if (score>0) and (score<high_threshold) and (score<best):
                    best = score
                    LGR.debug('Found better for {} in {}-{}'.format(
                        pdb, target_strand[j], target_strand[k+2]))
                    bestI['start1'] = target_strand[j]
                    bestI['start2'] = target_strand[j+2]
                    bestI['end1'] = target_strand[k]
                    bestI['end2'] = target_strand[k+2]
        LGR.info('Best:{}, {}-{}'.format(best, bestI['start1'], bestI['end2']))
        
        # Save the best alignment
        self.__pair_fit(bestI['start1'], bestI['start2'],
                        bestI['end1'], bestI['end2'])
        aligned = 'aligned'
        cmd.create(aligned, '{} and resi {}-{}'.format(self.targetN,
                                                       bestI['start1'],
                                                       bestI['end2']))
        name = '{}_{}.pdb'.format(prefix,
               os.path.basename(pdb).replace('.pdb','').replace('.gz', ''))
        cmd.save(os.path.join(self.save_path, name), aligned)

    def __call__(self):
        for pdb in self.pdbs:
            self.__find_best_match(pdb)


logging.basicConfig(level='DEBUG')
template = sys.argv[-2]
target = sys.argv[-1]
b = BestMatchFinder(template, target)
b()


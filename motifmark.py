#!/usr/bin/env python

##################################################################################################
## This program will take fasta file of genes (where exons are caps and introns are lower case) ##
## and a file of motifs (each motif on it's own line) and will produce a visualiztion of all    ##
## motif locations on each gene.                                                                ##
##################################################################################################

import cairo
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser("a program to produce two-line fastas")
    parser.add_argument("-f", "--fasta", type=str, help="fasta file of genes to pass through program", required=True)
    parser.add_argument("-m", "--motifs", type=str, help="motif file", required=True)
    parser.add_argument("-o", "--output", type=str, help="name of output (optional), default output uses name of fasta file", required=False)

    return parser.parse_args()
args = get_args()

###assign arguments to variables inside of program
fasta = args.fasta
motifs = args.motifs
out = args.output

#split input name for to get prefix to add to svg output
outputname = re.split('\.', fasta)

#motif disambiguation dictionary
regex_dict = {
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[TtUu]",
    "U":"[UuTt]",
    "W":"[AaTtUu]",
    "S":"[CcGg]",
    "M":"[AaCc]",
    "K":"[GgTtUu]",
    "R":"[AaGg]",
    "Y":"[CcTtUu]",
    "B":"[CcGgTtUu]",
    "D":"[AaGgTtUu]",
    "H":"[AaCcTtUu]",
    "V":"[AaCcGg]",
    "N":"[AaCcGgTtUu]",
    "Z":"[-]",
}

#############
# Functions #
#############

def draw_gene(sequence, gene_counter, header):
    '''   '''
    context.set_line_width(3)
    context.set_source_rgba(0, 0, 0, 1)
    context.move_to(15, int(gene_counter)*vert_pad + 30)
    context.line_to(15+len(sequence), int(gene_counter)*vert_pad + 30)
    context.stroke()
    context.move_to(15, int(gene_counter)*vert_pad + 10)
    context.show_text(header)
 

def find_exon(seq):
    '''this function will return the start and end positions of the exon'''
    exon_tuple = re.search('([A-Z]+)', seq)
    exon = exon_tuple.span()

    return exon

def draw_exon(exon, gene_counter):
    '''   '''
    context.set_source_rgba(0.23,0.25,0.25, 1)
    context.rectangle(15+int(exon[0]),int(gene_counter)*vert_pad+20,int(exon[1])-int(exon[0]),15)        #(x0,y0,x1,y1)
    context.fill()


def get_regex(motif):
    '''this function will take a motif and return a regex to find that motif'''
    motif_regex = ""
    for char in motif:
        motif_regex += regex_dict[char]

    return(motif_regex)

def find_positions(seq, reggie):
    '''   '''
    pos_list = []
    motifs = re.finditer(reggie, seq)
    for match in motifs:
        if match is None:
            continue
        else:
            coords = match.span()
            pos_list.append(coords)

    return pos_list

def draw_motif(R, G, B, coords, gene_counter):
    ''' '''
    context.set_source_rgba(float(R),float(G),float(B),.7)
    context.rectangle(15+int(coords[0]),int(gene_counter)*vert_pad+20,int(coords[1])-int(coords[0]),15)        #(x0,y0,x1,y1)
    context.fill()

#############
# Algorithm #
#############

#initialize gene count, empty string to build gene, and an empty list to store all of the genes
gene_count = 0
gene = ''
seqs = []

#iterate through the file to get the gene count and longest gene for svg dimensions
with open (fasta, "r") as fh:
    for line in fh:
        if line[0] != '>':
            seq = line.strip()
            gene += seq
        else:
            gene_count += 1
            seqs.append(gene)
            gene = ''
        seqs.append(gene)

longest_gene = (len(max(seqs, key=len)))

#svg set up using gene count:
vert_pad = 50 #this number dictates the space between genes, or the 'vertical padding'
width = int(longest_gene) + 30
height = int(gene_count+10) * vert_pad

#coordinates to display graphic and output name
surface = cairo.SVGSurface(outputname[0]+'.svg', width, height)
#create the coordinates you will be drawing on 
context = cairo.Context(surface)

motif_list = []
sequence = ''
gene_counter = 0
header_list = []

#establish motif colors
Rs = ("0.27","1","0.20","0.75","1")
Gs = ("0.5",".63","0.11","0.13",".72")
Bs = ("0.08",".15","1.00","0.06",".83")

with open (fasta, "r") as fh, open (motifs, "r") as mt:
#extract motifs from file into list
    for line in mt:
        motif = line.strip()
        motif = motif.upper()
        if motif not in motif_list:
            motif_list.append(motif)
#begin parsing genes for annotation
    for line in fh:
        if line[0] == '>':
            if sequence != '':
                draw_gene(sequence, gene_counter, header_list[gene_counter])
                exon = find_exon(sequence)
                draw_exon(exon, gene_counter)
                #make motif colors iterable
                itR = iter(Rs)
                itG = iter(Gs)
                itB = iter(Bs)
                for i in motif_list:
                    #call motif color
                    R = next(itR)
                    G = next(itG)
                    B = next(itB)
                    reggie = get_regex(i)
                    coords_list = find_positions(sequence, reggie)
                    for i in coords_list:
                        coords = list(i)
                        draw_motif(R,G,B,coords,gene_counter)
                gene_counter += 1
                sequence = ''
            header = line.strip()
            header_list.append(header)
        else:
            seq = line.strip()
            sequence += seq
    #repeat for last sequence
    draw_gene(sequence, gene_counter, header_list[gene_counter])
    exon = find_exon(sequence)
    draw_exon(exon, gene_counter)
    #make motif colors iterable again
    itR = iter(Rs)
    itG = iter(Gs)
    itB = iter(Bs)
    for i in motif_list:
        #call motif color
        R = next(itR)
        G = next(itG)
        B = next(itB)
        reggie = get_regex(i)
        coords_list = find_positions(sequence, reggie)
        for i in coords_list:
            coords = list(i)
            draw_motif(R,G,B,coords,gene_counter)
    #draw legend
    #make motif colors iterable one last time
    itR = iter(Rs)
    itG = iter(Gs)
    itB = iter(Bs)
    motif_counter = 0
    #write "legend"
    context.set_source_rgba(0, 0, 0, 1)
    context.move_to(15,int(gene_count)*vert_pad+15)
    context.show_text("Legend")
    for i in motif_list:
        #call motif color
        R = next(itR)
        G = next(itG)
        B = next(itB)
        #draw legend
        context.set_source_rgba(float(R), float(G), float(B), .9)
        context.rectangle(15,int(gene_count)*vert_pad+int(motif_counter+1)*20,40,10)
        context.fill()
        #write motif
        context.move_to(60,(int(gene_count)*vert_pad)+int((motif_counter+1)*20)+10)
        context.show_text(i)
        motif_counter += 1

surface.finish()
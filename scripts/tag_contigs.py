import os
import re

def cigar_tuple(cigar):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]
    tuples = []
    for i in range(len(lengths)):
        tuples.append([int(lengths[i]), ops[i]])
    return tuples

def main(snakemake):

    os.makedirs(f"results/{snakemake.wildcards.sample}/moddotplot/work/find_tigs/pafs/{snakemake.wildcards.hap}", exist_ok=True)
    os.makedirs(f"results/{snakemake.wildcards.sample}/moddotplot/work/find_tigs/beds/{snakemake.wildcards.hap}", exist_ok=True)
    centro_dict = {
            'chr13': {
                'p' : (15547593, 16522942),
                'q' : (16522942, 17498291),
            },
            'chr14': {
                'p' : (10092112, 11400261),
                'q' : (11400261, 12708411),
            },
            'chr15': {
                'p' : (16678794, 17186630),
                'q' : (17186630, 17694466),
            },
            'chr21': {
                'p' : (10962853, 11134529),
                'q' : (11134529, 11306205),
            },
            'chr22': {
                'p' : (12788180, 14249622),
                'q' : (14249622, 15711065),
            }
        }

    seq2len_dict = {}
    query2target2info_dict = {}
    qry_chrom_tracker = {}
    with open(snakemake.input.paf) as f:
        for line in f:
            try:
                query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]

                cg = [i.split(":")[-1] for i in line.strip().split("\t")[12:] if i[:2] == 'cg'][0]
                cg_tuple = cigar_tuple(cg)
                iden = round((sum([int(i[0]) for i in cg_tuple if i[1] == '=' or i[1] == 'M']) / sum(
                    [int(i[0]) for i in cg_tuple if i[1] in {'=', 'M', 'X', 'D', 'I'}])) * 100,2)

                query_len = int(query_len)
                query_start = int(query_start)
                query_end = int(query_end)
                target_len = int(target_len)
                target_start = int(target_start)
                target_end = int(target_end)

                seq2len_dict[query] = query_len
                seq2len_dict[target] = target_len

                if query not in qry_chrom_tracker:
                    qry_chrom_tracker[query] = []
                qry_chrom_tracker[query].append(target)

                # if int(alignment_len) < 200000:
                if int(alignment_len) <= 100000 or iden < 90:
                    continue
                if query not in query2target2info_dict:
                    query2target2info_dict[query] = {'p': [], 'q': []}

                if target_end - 1000000 >= centro_dict[target]['q'][1] and int(alignment_len) > 1000000:
                    arm = 'q'
                    query2target2info_dict[query][arm].append((
                    (query_start, query_end), strand, (target_start, target_end), target, num_matches, alignment_len, iden,
                    line.strip()))

                if target_start <= centro_dict[target]['p'][1]:
                    arm = 'p'
                    query2target2info_dict[query][arm].append((
                    (query_start, query_end), strand, (target_start, target_end), target, num_matches, alignment_len, iden,
                    line.strip()))
            except:
                continue


    fout_pq = open(snakemake.params.pq,'w')
    fout_p = open(snakemake.params.p, 'w')
    faln = open(snakemake.params.aln, 'w')
    for query, query2arm_dict in query2target2info_dict.items():
        try:
            query_len = seq2len_dict[query]
            ## contigs aligned to both p and q arm
            if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) > 0:
                distal_aln = sorted(query2arm_dict['q'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                qry_seq_start = 1
                qry_seq_end = distal_aln[0][1]
                target_chrom = distal_aln[3]
                if distal_aln[1] == '-':
                    qry_seq_start = distal_aln[0][0]
                    qry_seq_end = query_len
                # print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{query}',file=fout_pq)
                ## used for Arang's annotated assembly
                new_query = query.split('_')[1] + '_' + query.split('_')[0] if len(query.split('_')) > 1 else f'{query}_{target_chrom}'
                print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{snakemake.wildcards.sample}_{new_query}',file=fout_pq)
                # pq_contig_tracker.append(query)

            ## contig only aligned to p arm
            if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) == 0:
                distal_aln = sorted(query2arm_dict['p'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                qry_seq_start = 1
                qry_seq_end = distal_aln[0][1]
                target_chrom = distal_aln[3]
                if distal_aln[1] == '-':
                    qry_seq_start = distal_aln[0][0]
                    qry_seq_end = query_len
                # print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{snakemake.wildcards.sample}_{query}_{target_chrom}',file=fout_p)
                ## used for Arang's annotated assembly
                new_query = query.split('_')[1] + '_' + query.split('_')[0] if len(query.split('_')) > 1 else f'{query}_{target_chrom}'
                print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{snakemake.wildcards.sample}_{new_query}',file=fout_p)

            for arm, qry2arm_list in query2arm_dict.items():
                for aln_info in qry2arm_list:
                    print(aln_info[-1], file=faln)
        except:
            continue

    fout_pq.close()
    fout_p.close()
    faln.close()

    with open(str(snakemake.output.flag), "w") as f_flag:
        print(f"{snakemake.wildcards.sample}:{snakemake.wildcards.hap}:{snakemake.wildcards.region} Done.", file=f_flag)
if __name__=="__main__":
    main(snakemake)
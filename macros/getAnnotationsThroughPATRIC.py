import mackinac

mackinac.modelseed.ms_client.url = 'http://p3c.theseed.org/dev1/services/ProbModelSEED'
mackinac.workspace.ws_client.url = 'https://p3.theseed.org/services/Workspace'
mackinac.genome.patric_url = 'https://www.patricbrc.org/api/'
mackinac.get_token('asuccurro')


orgs = {'Root491': '1736548.3',
        'Root9': '1736604.3',
        'Root66D1': '1736582.3'}

mymodels = mackinac.list_modelseed_models()
mymodelsnames = ';'.join(map(lambda y: y['ref'], mymodels))


gprs = {'Root491': {},
        'Root9': {},
        'Root66D1': {}}

for o in ['Root491', 'Root9', 'Root66D1']:
#for o in ['Root491']:
    if orgs[o] not in mymodelsnames:
        print('Reconstructing modelseed model for', o, '(', orgs[o], ')')
        mackinac.get_genome_summary(orgs[o])
        mackinac.reconstruct_modelseed_model(orgs[o])
        mackinac.get_modelseed_model_stats(orgs[o])
    model = mackinac.get_modelseed_model_data(orgs[o])
    for r in model['modelreactions']:
        tmpgprs = []
        for p in r['modelReactionProteins']:
            for s in p['modelReactionProteinSubunits']:
                if s['feature_refs'] == []:
                    print(s)
                    gpr = 'Unknown'
                else:
                    gpr = ' OR '.join(map(lambda y: y.split('/')[-1], s['feature_refs']))
                tmpgprs.append('(%s)' % gpr)
        gprs[o][r['id']] = ' AND '.join(tmpgprs)

    with open('../data/MPI%s/MPI%s-gprs.tsv' % (o, o), 'w') as ofile:
        ofile.write('id\tgpr\n')
        for key, values in gprs[o].items():
            ofile.write('%s\t%s\n' % (key, values))
            
    with open('../data/MPI%s/MPI%s-reactions_with_GPR.tsv' % (o, o), 'w') as ofile:
        with open('../data/MPI%s/MPI%s-reactions.tsv' % (o, o), 'r') as ifile:
            for line in ifile:
                l = line.rstrip('\n')
                if 'gpr' in l:
                    ofile.write('%s\tgprfeat\n' % l)
                else:
                    rid = l.split('\t')[0]
                    ofile.write('%s\t%s\n' % (l, gprs[o].get(rid, 'Unknown')))


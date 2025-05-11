rename={}
with open('skin_project_meta.csv','r') as S:
    next(S)
    for i in S:
        i=i.split(',')[:4]
        rename[i[1]+'.sort.bam']=i[3]+i[2]+'.bam'

merge_ID={}
for i in rename:
    if not rename[i] in merge_ID:
        merge_ID[rename[i]]=[]
    merge_ID[rename[i]].append(i)

flag=0
for i in merge_ID:
    if len(merge_ID[i])==1:
        print('cp {} ../merged_bam/{} &'.format(merge_ID[i][0],i))
    else:
        if flag==4:
            print ('samtools merge ../merged_bam/{} {}'.format(i,' '.join(merge_ID[i])))
            flag=0
        else:
            print ('samtools merge ../merged_bam/{} {} &'.format(i,' '.join(merge_ID[i])))
            flag+=1


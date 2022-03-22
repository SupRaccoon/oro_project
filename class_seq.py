class SEQ:
    def __init__(self, data, iter_num):
        # название последовательности
        self.name = data['BlastOutput2']['report']['results']['search']['query_title']
        # получение результатов выравнивания
        self.ref = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][iter_num]['hseq']
        self.seq = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][iter_num]['qseq']
        self.qual = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][iter_num]['midline']
        # полная строка генома после выравнивания
        self.ref_full = str()
        self.seq_full = str()
        self.qual_full = str()
        # координаты старта и конца хота
        self.start = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][iter_num]['hit_from']
        self.end = data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'][iter_num]['hit_to']
        # словарь с нарезанными по координатам участками генома
        self.ref_seq_data = dict()
        self.new_seq_data = dict()
        self.qual_seq_data = dict()
        # списки мутаций
        self.snp_list = list()
        self.del_list = list()
        self.ins_list = list()
        # словарь с координатами CDS, словарь с количеством мутаций по каждому геноварианту, финальное заключение
        self.cord_cov = None
        self.mut_score = dict()
        self.curr_conclusion = list()
        # дополнительные переменные для определения сублиний отдельных геноварантов

    def complete_seq(self):

        # Функция для заполнения всего генома строками и пустыми полями, доводя до нужного размера
        if self.start > 1:
            self.ref_full = '-' * (self.start - 1) + self.ref
            self.seq_full = '-' * (self.start - 1) + self.seq
            self.qual_full = ' ' * (self.start - 1) + self.qual

        if self.end < 29903:
            self.ref_full = self.ref_full + '-' * (29903 - self.end)
            self.seq_full = self.seq_full + '-' * (29903 - self.end)
            self.qual_full = self.qual_full + ' ' * (29903 - self.end)

        self.seq_full = self.seq_full.upper()

    def get_cds_dict(self, cord_cov):
        # нарезка на кодирующие участки
        self.cord_cov = cord_cov
        for i in cord_cov.keys():
            start = cord_cov[i][0] - 1
            end = cord_cov[i][1]
            # start = start - self.start
            # end = end - self.start

            if set(self.seq_full[start:end]) != set('-'):

                self.ref_seq_data[i] = self.ref_full[start:end]
                self.new_seq_data[i] = self.seq_full[start:end]
                self.qual_seq_data[i] = self.qual_full[start:end]

    def find_cds(self, point):
        # сопоставление координаты и cds
        res = None
        for i in self.cord_cov.keys():
            start = self.cord_cov[i][0] - 1
            end = self.cord_cov[i][1]

            if start <= point <= end:
                res = i

                return res
        if not res:
            return str('none')

    def find_mutations(self):
        # заполнение переменных self.MUT_list
        if set(self.qual) != set('|'):
            id_mut = list()
            mut_pool = list()
            for i in range(len(self.qual)):
                if self.qual[i] != '|':
                    id_mut.append(i)

            for mut in id_mut:
                mut_pool.append([self.ref[mut], mut + self.start, self.seq[mut].upper(), self.find_cds(mut + self.start)])

            for elem in mut_pool:

                if elem[0] != '-' and elem[2] != '-':
                    self.snp_list.append(elem)
                elif elem[0] == '-' and elem[2] != '-':
                    self.ins_list.append(elem)

                elif elem[2] == '-' and elem[0] != '-':
                    self.del_list.append(elem)
                else:
                    print('smth wrong in finding mutations')

            # proc ins list
            ins_pool = list()
            if len(self.ins_list) > 1:
                doo = self.ins_list

                ins_start = None
                ins_end = None
                nuc = None
                for id_ins in range(len(doo)-1):
                    if doo[id_ins][1] + 1 == doo[id_ins + 1][1]:
                        if ins_start:
                            ins_end += 1
                            nuc += doo[id_ins + 1][2]
                        else:
                            ins_start = doo[id_ins][1]
                            nuc = doo[id_ins][2] + doo[id_ins + 1][2]
                            ins_end = doo[id_ins + 1][1]
                    else:
                        if ins_start:
                            ins_pool.append(['ins', ins_start, nuc, doo[id_ins][3]])
                        else:

                            ins_pool.append(['ins', doo[id_ins][1], doo[id_ins][2], doo[id_ins][3]])

                        ins_start = None
                        ins_end = None
                        nuc = None
                if ins_start:
                    ins_pool.append(['ins', ins_start, nuc, doo[-1][3]])

            # proc ins coord
            if len(self.ins_list) > 0:

                for ins in self.ins_list:
                    for snp in self.snp_list:
                        if ins[1] < snp[1]:
                            snp[1] -= 1
                    for deletion in self.del_list:
                        if ins[1] < deletion[1]:
                            deletion[1] -= 1

            if len(ins_pool) > 0:
                for ins in ins_pool:

                    for insertion in ins_pool:
                        if ins[1] < insertion[1]:
                            insertion[1] -= 1
                self.ins_list = ins_pool

            elif len(self.ins_list) == 1:
                self.ins_list = [['ins', self.ins_list[0][1], self.ins_list[0][2], self.ins_list[0][3]]]
            # proc del list

            if len(self.del_list) > 1:
                foo = self.del_list

                del_pool = list()
                del_start = None
                del_end = None
                nuc = None
                for id_del in range(len(foo)-1):
                    if foo[id_del][1] + 1 == foo[id_del+1][1]:
                        if del_start:
                            del_end += 1
                            nuc += foo[id_del+1][0]
                        else:
                            del_start = foo[id_del][1]
                            nuc = foo[id_del][0] + foo[id_del+1][0]
                            del_end = foo[id_del+1][1]
                    else:
                        if del_start:
                            del_pool.append(['del', '{} - {}'.format(del_start, del_end), nuc, foo[id_del][3]])
                        else:

                            del_pool.append(['del', str(foo[id_del][1]), foo[id_del][0], foo[id_del][3]])

                        del_start = None
                        del_end = None
                        nuc = None
                del_pool.append(['del', '{} - {}'.format(del_start, del_end), nuc, foo[-1][3]])

                self.del_list = del_pool
        #
        # else:
        #     print('zero mutations found')

    def compare_mut(self, curr_mut_list):
        # поиск конкретных мутаций каждой клады, на выходе массив со словарями формата
        # {'mut': {'L5F': 'missing', 'T95I': 'missing', 'D253G': 'missing'}, 'count': [0, 3], 'label': 'Iota'}
        DEL = curr_mut_list['del']
        INS = curr_mut_list['ins']
        SNP = curr_mut_list['snp']

        max_count = curr_mut_list['count']
        label = curr_mut_list['label']
        found = list()
        not_found = list()
        un_cover = list()

        for elem in DEL:
            if elem[0] in self.del_list:
                found.append(elem[1])
            elif self.start <= int(elem[0][1].split(' - ')[0]) <= self.end:
                not_found.append(elem[1])
            else:
                un_cover.append(elem[1])
        for elem in INS:

            if elem[0] in self.ins_list:
                found.append(elem[1])
            elif self.start <= elem[0][1] <= self.end:
                not_found.append(elem[1])
            else:
                un_cover.append(elem[1])
        for elem in SNP:
            if elem[0] in self.snp_list:
                found.append(elem[1])
            elif self.start <= elem[0][1] <= self.end:
                if elem[0][2] != 'N':
                    not_found.append(elem[1])
                else:
                    un_cover.append(elem[1])
            else:
                un_cover.append(elem[1])
        mut = dict()
        for elem in found:
            mut[elem] = 'found'
        for elem in not_found:
            mut[elem] = 'not found'
        for elem in un_cover:
            mut[elem] = 'missing'

        report = dict()
        report['mut'] = mut
        report['count'] = [len(found), max_count]
        report['label'] = label
        # print(report['mut'])
        # print()
        return report

    def unite_mut(self, report):
        # да много pass но так удобнее ориентироваться

        # на вход массив из функции compare_mut и сам объект (один хит, то есть либо часть сиквенса, либо весь),
        # результат : большой словарь по всем кладам

        self.mut_score['mut'] = dict()
        self.mut_score['label'] = dict()
        for rep in report:
            m_pool = rep['mut']
            self.mut_score['label'][rep['label']] = rep['count']

            for mut in m_pool.keys():
                if m_pool[mut] == 'found':
                    self.mut_score['mut'][mut] = 'found'
                elif mut in self.mut_score['mut'].keys():
                    if self.mut_score['mut'][mut] == 'found':
                        pass
                    elif self.mut_score['mut'][mut] == 'not found':
                        pass
                    elif self.mut_score['mut'][mut] == 'missing':
                        if m_pool[mut] == 'not found':
                            self.mut_score['mut'][mut] = m_pool[mut]
                        else:
                            pass
                else:
                    self.mut_score['mut'][mut] = m_pool[mut]

    def unite_seq(self, another_seq, list_of_mutations):
        # да много pass но так удобнее ориентироваться

        # объединяет self.mut_score для нескольких объектов
        m_pool = another_seq.mut_score['mut']

        copy_pool = []
        for mut in m_pool.keys():
            if m_pool[mut] == 'found':

                if self.mut_score['mut'][mut] == 'found':
                    copy_pool.append(mut)
                self.mut_score['mut'][mut] = 'found'
            elif mut in self.mut_score['mut'].keys():
                if self.mut_score['mut'][mut] == 'found':
                    pass
                elif self.mut_score['mut'][mut] == 'not found':
                    pass
                elif self.mut_score['mut'][mut] == 'missing':
                    if m_pool[mut] == 'not found':
                        self.mut_score['mut'][mut] = m_pool[mut]
                    else:
                        pass
            else:
                self.mut_score['mut'][mut] = m_pool[mut]

        an_count = another_seq.mut_score['label']
        for elem in an_count.keys():
            self.mut_score['label'][elem][0] += an_count[elem][0]

        for elem in copy_pool:
            for clade in list_of_mutations.keys():
                if elem in list_of_mutations[clade]:
                    self.mut_score['label'][clade][0] -= 1
        self.snp_list += another_seq.snp_list
        self.del_list += another_seq.del_list
        self.ins_list += another_seq.ins_list

    def conclusion_clade(self, list_of_mutations, count_list):

        set_uncoverage = set()
        for clade in list_of_mutations.keys():
            all_count = count_list[clade]

            f_count = 0
            nf_count = 0
            nc_count = 0
            # g_count = 0

            foo = False
            # print(self.mut_score)
            for mut in list_of_mutations[clade]:
                if self.mut_score['mut'][mut] == 'found':
                    f_count += 1
                    foo = True
                    # print('found_mut', self.mut_score['mut'][mut])
                elif self.mut_score['mut'][mut] == 'missing':
                    nc_count += 1
                else:
                    nf_count += 1
            N = 6
            g_count = count_list[clade] - nc_count
            gap_for_conclusion = 1 + (g_count - 3)//N

            # print(g_count, nc_count, len(list_of_mutations[clade]), all_count, 'before')
            # print(clade, f_count, nf_count, nc_count)
            if 'недостаточное покрытие' not in self.curr_conclusion:
                if f_count > 2:
                    print(clade, 'found:', f_count, 'not found', nf_count, 'gap', gap_for_conclusion)
                    if nf_count < gap_for_conclusion:
                        self.curr_conclusion.append(clade)
                    elif nf_count < 2 * gap_for_conclusion:
                        self.curr_conclusion.append('пред. ' + clade)
                elif f_count == 2:
                    if nf_count <= 1 and foo:
                        self.curr_conclusion.append('пред. ' + clade)
                elif all_count - nc_count < 2:
                    set_uncoverage.add(clade)
            # print(g_count, nc_count, len(list_of_mutations[clade]), clade)
        if len(self.curr_conclusion) == 0:
            for voc in ["Alpha", "Beta", "Gamma", "Delta"]:
                if voc in set_uncoverage:
                    self.curr_conclusion = ['недостаточное покрытие']
                    # print(clade)
                    #     print(g_count, nc_count, len(list_of_mutations[clade]))
            # print(set_uncoverage)

        if len(self.curr_conclusion) == 0:
            self.curr_conclusion.append('иной')

    def find_lineage(self, lineage_pool):
        # по идее это должно выглядеть так:
        # на вход уже собранный объект и если он относится к группе с доп линиями,
        # то проходит определение конкретной линии
        # на выходе будет доп атрибут который меняет "заключение" на конкретную линию или комбинацию если
        # такие были найдены
        # TODO: 1. добавить json с мутациями для каждой линии
        #       2. реализовать поиск целевых мутаций и выставление нового заключения
        #       3. найти штук 5-6 геномов для тестов и протестировать
        pass

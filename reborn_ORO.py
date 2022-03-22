import os
import uuid
import pandas as pd
from tqdm.auto import tqdm  # new lib
import inspect
import argparse
import json
# import sys
import tempfile
from class_seq import SEQ

# TODO: 1. Add change way to start blast
#       2. Add Delta lineages


def make_table(curr_list_of_seq, curr_config):
    """
    Создание уже вспомогательной таблицы для отчёта
    :param curr_list_of_seq: массив объектов класса SEQ после анализа
    :param curr_config: конфиги из папки data
    :return: таблица формата "имя"/"все мутации"/"заключение"
    """
    df_d = pd.DataFrame(columns=['Name'] + curr_config['main_table'] + ['Conclusion'])
    for sequence in curr_list_of_seq:
        pre_table = [sequence.name]
        for _mut in curr_config['main_table']:
            pre_table.append(sequence.mut_score['mut'][_mut])

        concl = sequence.curr_conclusion[0]
        for elem in sequence.curr_conclusion[1:]:
            concl += ' или ' + elem
        pre_table.append(concl)
        df_d.loc[len(df_d)] = pre_table

    return df_d


def to_str(ls):
    """
    Вспомогательная функция для перевода массива в строку
    :param ls: массив ["мутация_1", ..., "мутация_n"]
    :return: строка "мутация_1, ..., мутация_n"
    """
    c = ''
    if len(ls) > 0:
        for _c in range(len(ls) - 1):
            c += ls[_c] + ', '
        c += ls[-1]
    return c


def make_clade_table(curr_list_of_seq, curr_config, curr_gene_variant):
    """
    По сути то же самое что и "make_table", только по каждому геноварианту
    :param curr_list_of_seq: массив объектов класса SEQ после анализа
    :param curr_config: конфиги из папки data
    :param curr_gene_variant: геновариант для которого делается таблица
    :return:
    """
    df_d = pd.DataFrame(columns=['Name'] + curr_config['clade'][curr_gene_variant] + ['Conclusion'])
    for sequence in curr_list_of_seq:
        pre_table = [sequence.name]
        for _mut in curr_config['clade'][curr_gene_variant]:
            pre_table.append(sequence.mut_score['mut'][_mut])

        concl = sequence.curr_conclusion[0]
        for elem in sequence.curr_conclusion[1:]:
            concl += ' или ' + elem
        if curr_gene_variant in concl.split(' '):
            pre_table.append(concl)
        else:
            pre_table.append('another variant')
        df_d.loc[len(df_d)] = pre_table

    return df_d


def make_all_table(curr_list_of_seq, cov_list):
    """
    Создание полной таблицы с мутациями
    :param curr_list_of_seq: массив объектов класса SEQ после анализа
    :param cov_list: данные по покрытию для каждого образца
    :return: таблица со всеми мутациями
    """
    df_d = pd.DataFrame(columns=['ref', 'loc', 'seq', 'cds'])
    for sequence in curr_list_of_seq:
        pre_table = [sequence.name, 0, '', '']
        df_d.loc[len(df_d)] = pre_table

        for elem in range(len(sequence.snp_list)):
            sequence.snp_list[elem] = [sequence.snp_list[elem][0], str(sequence.snp_list[elem][1]),
                                       sequence.snp_list[elem][2], sequence.snp_list[elem][3]]
        for elem in range(len(sequence.del_list)):
            if sequence.del_list[elem][1] != 'None - None':
                sequence.del_list[elem] = [sequence.del_list[elem][0], str(sequence.del_list[elem][1]),
                                           sequence.del_list[elem][2], sequence.del_list[elem][3]]
        for elem in range(len(sequence.ins_list)):

            sequence.ins_list[elem] = [sequence.ins_list[elem][0], str(sequence.ins_list[elem][1]),
                                       sequence.ins_list[elem][2], sequence.ins_list[elem][3]]

        MUT = []

        for elem in sequence.snp_list:
            if elem not in MUT and elem[2] != 'N':
                MUT.append(elem)
        for elem in sequence.del_list:
            if elem not in MUT and elem[1] != 'None - None':
                MUT.append(elem)
        for elem in sequence.ins_list:
            if elem not in MUT:
                MUT.append(elem)

        for j in range(len(MUT) - 1):
            for ii in range(len(MUT) - j - 1):
                if int(MUT[ii][1].split('-')[0]) > int(MUT[ii + 1][1].split('-')[0]):
                    MUT[ii], MUT[ii + 1] = MUT[ii + 1], MUT[ii]

        for elem in MUT:
            df_d.loc[len(df_d)] = elem
    cov_str = ''
    for gap in cov_list:
        cov_str += f'{gap[0]} - {gap[1]},'
    df_d.loc[len(df_d)] = ['No coverage', cov_str, '', '']
    return df_d


def check_coverage(gap_pool):
    """
    Функция для создания списка непокрытых областей
    :param gap_pool: массив покрытия из бласта
    :return: массив непокрытых областей
    """
    pool = sorted(gap_pool, key=lambda elem: elem[0])
    curr_cov = [pool[0]]
    for i in range(len(pool)):
        if pool[i][0] > curr_cov[-1][1]:
            curr_cov.append(pool[i])
        elif pool[i][1] > curr_cov[-1][1]:
            curr_cov[-1][1] = pool[i][1]

    result = [[1, curr_cov[0][0] - 1]]
    for i in range(len(curr_cov)-1):
        result.append([curr_cov[i][1] + 1, curr_cov[i+1][0] - 1])
    result.append([curr_cov[-1][1] + 1, 29903])

    return result


def parus(args_data, args_result, args_one_sample, args_table=False, args_gapopen=1, args_gapextend=1):
    """
    Осноная функция скрипта
    :param args_data:       путь к файлу с последовательностью (или несколькими в одном файле)
    :param args_result:     префикс для результирующих файлов
    :param args_one_sample: флаг, передаётся в случае если все последовательности в файле относятся к одному образцу
                            то есть, если там один полный геном то от передачи этого флага ничего не поменяется, а вот
                            фрагменты будут анализироваться либо как отдельные последовательности, либо как части одной
    :param args_table:      флаг, если передан, то се таблицы сохраняются, если нет, то после сохранения в финальный
                            словарь удаляются
    :param args_gapopen:    параметры для бласта (указаны по умолчанию, и можно забить)
    :param args_gapextend:  параметры для бласта (указаны по умолчанию, и можно забить)
    :return:                словарь с примерно такой структурой: ключи -- ("result", "raw_all_table", "raw_report_table"),
                            первый -- словарь с результатами анализа,
                            второй -- csv со всеми мутациями строкой,
                            третий -- csv с отчётной таблицей строкой
    """

    # определение в какой папке лежит скрипт
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    path_file = os.path.dirname(os.path.abspath(filename))

    # получение путей ко всем необходимым файлам (конфиги, данные, и база для бласта)
    path_to_ref_seq = os.path.join(path_file, 'data/ref_seq.fasta')
    path_to_curr_file = args_data
    path_to_data = os.path.join(path_file, 'data')
    path_to_list_of_clade = os.path.join(path_file, 'clade')
    list_of_clade = os.listdir(path_to_list_of_clade)

    # запуск бласта
    label_result_file = args_result
    label_tmp_file = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
    blast_path = os.path.dirname(os.path.realpath(__file__)) + '/tools/blastn'
    os.system('chmod +x %s' % blast_path)
    run_blast = '{} -db {} -query {} -out {} -gapopen {} -gapextend {} -outfmt 13'.format(blast_path, path_to_ref_seq,
                                                                                          path_to_curr_file, label_tmp_file, args_gapopen, args_gapextend)

    os.system(run_blast)
    # os.system(
    #     'blastn -db {}  -query {}  -out {} -gapopen 1 -gapextend 1 -outfmt 13'.format(path_to_ref_seq,
    #                                                                                   path_to_curr_file,
    #                                                                                   label_tmp_file))

    # чтение конфига
    with open(os.path.join(path_to_data, 'config.json'), 'r') as f:
        config = json.load(f)
    
    # чтение результатов бласта и обработка ошибок
    try:
        with open(label_tmp_file, 'r') as f:
            result = json.load(f)
    except:
        rep_file = dict()
        rep_file["result"] = [dict()]
        rep_file["result"][0]['conclusion'] = ["ERROR"]
        rep_file["result"][0]['found'] = ['-']
        rep_file["result"][0]['not_found'] = ['-']
        rep_file["result"][0]['missing'] = ['-']
        rep_file["result"][0]['name'] = ["ERROR"]
        rep_file["raw_all_table"] = ''
        rep_file["raw_report_table"] = ''

        return rep_file

    # если всё ок сборка результатов бласта в одно место
    pool_seq = list()
    for path in result['BlastJSON']:
        with open(os.path.join(tempfile.gettempdir(), path['File']), 'r') as f:
            pool_seq.append(json.load(f))

    # чтение ещё одного сопутствующего словаря с координатами кодирующих частей генома
    with open(os.path.join(path_to_data, 'cord_cov.json'), 'r') as f:
        cord_cov = json.load(f)

    multiple_seq = list()
    list_of_seq = list()
    start_end_list = list()

    # начало самого анализ
    for data in tqdm(pool_seq):
        one_seq = list()
        # очередная обработка ошибок
        if 'hits' in data['BlastOutput2']['report']['results']['search']:
            if len(data['BlastOutput2']['report']['results']['search']['hits']) == 0:
                rep_file = dict()
                rep_file["result"] = [dict()]
                rep_file["result"][0]['conclusion'] = ["ERROR"]
                rep_file["result"][0]['found'] = ['-']
                rep_file["result"][0]['not_found'] = ['-']
                rep_file["result"][0]['missing'] = ['-']
                rep_file["result"][0]['name'] = ["ERROR"]
                rep_file["raw_all_table"] = ''
                rep_file["raw_report_table"] = ''

                return rep_file
            else:
                pass
        else:
            rep_file = dict()
            rep_file["result"] = [dict()]
            rep_file["result"][0]['conclusion'] = ["ERROR"]
            rep_file["result"][0]['found'] = ['-']
            rep_file["result"][0]['not_found'] = ['-']
            rep_file["result"][0]['missing'] = ['-']
            rep_file["result"][0]['name'] = ["ERROR"]
            rep_file["raw_all_table"] = ''
            rep_file["raw_report_table"] = ''

            return rep_file

        for curr_iter in range(len(data['BlastOutput2']['report']['results']['search']['hits'][0]['hsps'])):
    
            seq = SEQ(data, curr_iter)
            start_end_list.append([seq.start, seq.end])
            seq.complete_seq()
            seq.get_cds_dict(cord_cov)
            seq.find_mutations()
    
            pool_of_clade = list()

            # чтение доп файлов с мутациями по каждому из геновариантов
            for clade in list_of_clade:
                with open(os.path.join(path_to_list_of_clade, clade), 'r') as f:
                    curr_clade = json.load(f)
                pool_of_clade.append(seq.compare_mut(curr_clade))
    
            seq.unite_mut(pool_of_clade)
            one_seq.append(seq)

        if len(one_seq) > 1:
            for i in range(1, len(one_seq)):
                one_seq[0].unite_seq(one_seq[i], config['clade'])
        if not args_one_sample:
            one_seq[0].conclusion_clade(config['clade'], config['count'])
        multiple_seq.append(one_seq[0])

    if args_one_sample:
        if len(multiple_seq) > 1:
            for i in range(1, len(multiple_seq)):
                multiple_seq[0].unite_seq(multiple_seq[i], config['clade'])
        multiple_seq[0].conclusion_clade(config['clade'], config['count'])
    
        list_of_seq.append(multiple_seq[0])
    else:
        list_of_seq = multiple_seq

    # удаление результатов бласта после анализа
    os.system('rm {}*'.format(label_tmp_file))
    # сборка финальных таблиц
    df = make_table(list_of_seq, config)
    all_df = make_all_table(list_of_seq, check_coverage(start_end_list))

    new_df = pd.DataFrame(columns=["Name", "Found", "Not found", "No coverage", "Conclusion"])

    file_for_vgarus = dict()
    file_for_vgarus['result'] = []

    for name in df.Name.values:
        file_for_vgarus['result'].append(dict())
        good = []
        good_js = []
        not_good = []
        not_good_js = []
        idk = []
        idk_js = []
        row = df[df.Name == name]

        for mut in row.columns[1:-1]:

            if (row[mut] == "found").values[0]:
                good.append(config['label_for_table'][mut])
                good_js.append(mut)

            elif (row[mut] == "not found").values[0]:
                not_good.append(config['label_for_table'][mut])
                not_good_js.append(mut)
            else:
                idk.append(config['label_for_table'][mut])
                idk_js.append(mut)
        print(row.Conclusion.values[0])
        file_for_vgarus['result'][-1]['conclusion'] = row.Conclusion.values[0]
        file_for_vgarus['result'][-1]['found'] = good_js
        file_for_vgarus['result'][-1]['not_found'] = not_good_js
        file_for_vgarus['result'][-1]['missing'] = idk_js
        file_for_vgarus['result'][-1]['name'] = name
        curr_row = [name, to_str(good), to_str(not_good), to_str(idk), row.Conclusion.values[0]]
        new_df.loc[len(new_df)] = curr_row

    # сохраняем чтобы потом перевести в строку
    all_df.to_csv(os.path.join(tempfile.gettempdir(), label_result_file + '_all' + '.csv'), index=False)
    new_df.T.to_csv(os.path.join(tempfile.gettempdir(), label_result_file + '_report' + '.csv'), header=False)

    # собираем всё вместе в один словарик
    if args_table:
        for clade in config['clade'].keys():
            make_clade_table(list_of_seq, config, clade).to_csv(label_result_file + '_' + clade + '.csv', index=False)

    else:
        with open(os.path.join(tempfile.gettempdir(), label_result_file + '_all' + '.csv'), 'r') as f:
            file_for_vgarus['raw_all_table'] = f.read()
        with open(os.path.join(tempfile.gettempdir(), label_result_file + '_report' + '.csv'), 'r') as fr:
            file_for_vgarus['raw_report_table'] = fr.read()

        os.system(f" rm {os.path.join(tempfile.gettempdir(), label_result_file + '_all' + '.csv')} ")
        os.system(f" rm {os.path.join(tempfile.gettempdir(), label_result_file + '_report' + '.csv')} ")
        return file_for_vgarus


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data", help="full path to file", type=str)
    parser.add_argument("-result", "--result", help="Name of the result file", type=str, default='result')
    # parser.add_argument('-no_code', '--no_code', action='store_true', help="Change standard unicode (UTF-8) to cyrillic encode", required=False)
    parser.add_argument('-one_sample', '--one_sample', action='store_true', help="Use pooling", required=False)
    parser.add_argument('-table', '--table', action='store_true', help="Add full pack of table",
                        required=False)
    parser.add_argument('-gapopen', '--gapopen',  help="gapopen blust",
                        required=False, default=1)
    parser.add_argument('-gapextend', '--gapextend',  help="gapextend blust",
                        required=False, default=1)
    args = parser.parse_args()

    print(parus(args.data, args.result, args.one_sample, args.table, int(args.gapopen), int(args.gapextend)))

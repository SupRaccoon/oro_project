import json
import os
from threading import Thread
from queue import Queue
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
try:
    from tqdm.auto import tqdm
except ModuleNotFoundError:
    print("As the additional script was run, please, install everything from the requirements.txt, " +
          "and then rerun this script. Not all the modules were installed via setup when the 'faquc' module " +
          "was installed, but only the required for the 'faquc'.")
    exit()

from reborn_ORO import parus


def worker(path_data, fasta_name):
    """
    Функция для выполнения каждым из потоков ThreadPool. Считывает fasta-файл, запускает пайплайн, возвращает
    результат работы пайплайна. В функции предусмотрена заглушка на случай получения пустого имени файла. \n \n
    """
    # if len(fasta_name) > 0:
    #     # считываем файл, чтобы в пайплайн передать fasta-like строку
    #     with open(fasta_name, "r", encoding="utf-8") as fao:
    #         ff_str = fao.read()
    # else:
    #     ff_str = ""
    # конфигуратор используется 'глобально', поэтому можно сюда передавать без ключей
    # result = faquco.favalid.pipeline.run_pipeline(ff_str, config['pipeline'], ff_id)

    return parus(os.path.join(path_data, fasta_name), f'{fasta_name.split(".")[0]}', True, False)


def consume():
    """
    Функция, выполняемая потоком записи результатов. Проверяет наличие результатов работы скрипта в очереди,
    пишет в файлы. Запись идёт в два файла - один файл TSV для удобного импорта через DataFrame,
    другой файл хранит результат работы пайплайна в сыром виде.

    Заканчивает работу в случае, если в основном теле скрипта кончилась обработка файлов.
    """
    while True:
        if not queue.empty():
            (sub_name, body) = queue.get()
            # пишем в обработанном виде только нужные поля для считывания как tsv
            # result_out_tsv.write(f"{sub_name}\t{body.get('file_valid')}\t{body.get('code_style')}" +
            #                      f"\t{body.get('quality_group')}\t{str(body.get('reasons'))}\n")
            # # пишем весь результат работы в сыром виде в txt файл
            result_out_txt.write(f"{sub_name}\t{json.dumps(body, ensure_ascii=False)}\n")
        if ind == NUM_OF_FILES - 1:
            # result_out_tsv.close()
            result_out_txt.close()
            return


def errors_consumer():
    """
    Функция, выполняемая потоком записи ошибок. Проверяет наличие ошибок работы скрипта в очереди, пишет в файл
    с так сказать логами.

    Заканчивает работу в случае, если в основном теле скрипта кончилась обработка файлов.
    """
    while True:
        if not queue_2.empty():
            message = queue_2.get()
            errors_out_txt.write(message + "\n")  # пишем весь результат работы в сыром виде в txt файл
        if ind == NUM_OF_FILES - 1:
            errors_out_txt.close()
            return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data", help="full path to file", type=str)
    parser.add_argument("-result", "--result", help="Name of the result file", type=str, default='result')
    # parser.add_argument('-no_code', '--no_code', action='store_true', help="Change standard unicode (UTF-8) to cyrillic encode", required=False)

    args = parser.parse_args()  # парсим аргументы

    # получаем непустую директорию с fasta-файлами
    try:
        PATH_FASTA = os.listdir(args.input)
        NUM_OF_FILES = len(PATH_FASTA)
        assert len(PATH_FASTA) != 0
    except FileNotFoundError:
        print("Input directory was not found")
        exit()
    except AssertionError:
        print("Input directory is empty")
        exit()
    else:
        # если все ок, то продолжаем
        ind = 0  # объявляем переменную номера обрабатываемого файла перед созданием очередей

        # очередь для записи результатов работы пайплайна
        queue = Queue()
        # result_out_tsv = open(args.out_name + ".tsv", "w", encoding="utf-8")  # тут используется имя из ключей
        result_out_txt = open(args.out_name + ".txt", "w", encoding="utf-8")  # тут используется имя из ключей
        # очередь для записи ошибок скрипта
        queue_2 = Queue()
        errors_out_txt = open(args.out_name + "errors.txt", "w", encoding="utf-8")
        # поток для обработки первой очереди
        consumer = Thread(target=consume)
        consumer.setDaemon(True)
        consumer.start()
        # поток для обработки второй очереди
        err_writer = Thread(target=errors_consumer)
        err_writer.setDaemon(True)
        err_writer.start()

        # config = faquco.config_once()
        # запускаем ThreadPool
        with ThreadPoolExecutor(max_workers=args.num_threads) as executor:
            # добавляем progress bar для оценки скорости
            progress_bar = tqdm(total=NUM_OF_FILES)
            # создаем словарь объектов запуска
            future_to_file = {executor.submit(worker,
                                              args.input + name,
                                              id_ed): name for id_ed, name in enumerate(PATH_FASTA)}
            # итерируемся (с индексом) по полученному словарю
            for ind, future in enumerate(as_completed(future_to_file)):
                f_name = future_to_file[future]
                try:
                    # если проблем нет, то просто получаем результат
                    data = future.result()
                except Exception as exc:
                    # если проблемы есть, то передаем их в очередь для записи ошибок
                    queue_2.put('%r generated an exception: %s' % (f_name, exc))
                else:
                    # который передаем в очередь для результатов
                    queue.put((os.path.splitext(f_name)[0], data))
                progress_bar.update()  # на каждый файл обновляем статус-бар
            progress_bar.close()
        # ожидаем завершения работы потоков, работающих с очередями
        consumer.join()
        err_writer.join()
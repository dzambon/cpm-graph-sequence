for ds in ['dog{}'.format(i) for i in range(1, 5)] + ['hum{}'.format(i) for i in range(1, 9)]: 
    with open('datasets/kaggle-seizure/{}/info.txt'.format(ds)) as f:
        print(f.read())
    f.close()

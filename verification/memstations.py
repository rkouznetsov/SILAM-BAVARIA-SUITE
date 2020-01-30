
# -*- coding: utf-8 -*-
from __future__ import print_function
from toolbox import timeseries

import pytz


ignore=[u'М1 ', u'М2 ', u'М3 ', u'М4 ', u'М5 ', u'М6 ',]

name2code={
    u'Балчуг':'Balchug',
    u'Бирюлево':'Birulevo',
    u'Бутлерова':'Butlerova',
    u'Чаянова':'Chayanova',
    u'Черемушки':'Cheremushki',
    u'Долгопрудная':'Dolgoprudnaya',
    u'Глебовская':'Glebovskaya',
    u'Гурьевский проезд':'Gurievskiy',
    u'Гурьянова':'Guriyanova',
    u'Хамовники':'Hamovniki',
    u'Казакова':'Kazakova',
    u'Новокосино':'Kosino',
    u'Кожухово':'Kozhuhovo',
    u'Лосиный остров':'Losiniy_ostr',
    u'Люблино':'Lublino',
    u'МАДИ':'Madi',
    u'Марьино':'Marino',
    u'Головачева':'MNPZGolovacheva',
    u'Капотня':'MNPZKapotnya',
    u'МГУ':'MSU',
    u'Народного ополчения':'Narodnogo_Opolcheniya',
    u'Нижняя Масловка':'Nizhnyaya_Maslovka',
    u'МАДИ':'Madi',
    u'Останкино 0':'Ostankino_0',
    u'Площадь Гагарина':'Ploshad_Gagarina',
    u'Полярная ул.':'Polyarnaya_ul',
    u'пос.Кузнецово':'Pos_Kuznetsovo',
    u'Пролетарский проспект':'Proletarskiy_Prospect',
    u'Рогово':'Rogovo',
    u'Семенково':'Semenkovo',
    u'Шаболовка':'Shabolovka',
    u'Спартаковская пл.':'Spartakovskaya',
    u'Спиридоновка ул.':'Spiridonovka',
    u'Сухаревская пл.':'Suharevskaya_pl',
    u'Толбухина':'Tolbuhina',
    u'Троицк':'Troick',
    u'Вешняки':'Veshnyaki',
    u'Зеленоград 11':'Zelenograd_11',
    u'Зеленоград 16':'Zelenograd_16',
    u'Зеленоград 6':'Zelenograd_6',
    u'Туристская':'Turistskaya',
    u'Ак. Анохина':'Ak_anokhina',
    u'Звенигород':'Zvenigorod',
    u'Кожуховский проезд':'Kozhuhovskiy_proezd',
    u'Коптевский бул.':'Koptevskiy_bul',
    u'МКАД 105 км (восток)':'MKAD_105',
    u'Мелитопольская':'Melitopolskaya'}

localtz=pytz.timezone('Europe/Moscow')

def getStations(tsfile):

    print ("Getting stations", tsfile)
    stations=[]
    with open(tsfile) as inf:
        for l in inf:
            a=l.split()
    #        Balchug 55.7453 37.627975 127  Roadside
            
            stations.append(timeseries.Station(a[0], a[0], a[2], a[1], a[3], a[4]))
    #        code, name, x, y, height=0.0, area_type='', dominant_source='')
    stationdic={}
    for st in stations:
        stationdic[st.code] = st
    
    return stationdic

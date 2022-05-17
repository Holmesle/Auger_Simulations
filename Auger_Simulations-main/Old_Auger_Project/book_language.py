import numpy as np

pawnee_pre = ['wist','bit','hiw','diss','diit','shuuw','sak','niish','kuuk','hak','haht',\
            'yaas','baws','wkib','nayaw','kah',"ka'an",'kab', 'had', "hi'uus",'ban','dahk', 'diit'\
            'tan', 'tah',"n'aw", 'ut','wih',"ts'at",'kik','bat',"cha'y",'kid']

pawnee_suf = ['itsi','amah','aha',"a'",'isi','ihit','iwah','at','ahwin','ush','akuh'\
            'ayah','idu','ikah','ayu','usuh','iku', 'inu', 'itsi','iwi','kian','attih'\
            'asa','awsa']


# prefix = np.random.choice(pawnee_pre)
# sufix = np.random.choice(pawnee_suf)

# print(prefix + sufix)

Uk_pre = ['birm','sheff','camb','manch','glas','edin','dund'\
    'aber', 'liver','dub','lim','brigh','lond','east','south'\
    'brit','card','nott','peter','chest','black','james','rob'\
    'will','rich','nor','abin','alf','als','ash','ban','barns'\
    'beac','bed','berk','bing','bridge','buck','castle','ches'\
    'cinder', 'dart','dews','don','dun','durs','farn','finch',\
    'gray','hales','harl','kimb','kirk','leigh','loft','long',\
    'mad','marl','nel','pain','pet','pres','ray','roch','roth'\
    'silver','strat','swan','thom','thorp','tod','uck','ox',\
    'ver','wall','war','wat']

Uk_suf = ['ington', 'asey', 'ingham', 'mouth', 'bridge', 'irsk', 'caster'
          'worth', 'nely', 'ton', 'inster', 'ford', 'anley', 'borough', 'land'
          'en', 'bourne', 'burry', 'brook', 'ton', 'elsey', 'field', 'bergh'
          'down', 'mond', 'hurst', 'bach', 'combe', 'ston', 'well', 'chester', 'barns'
          'dale', 'gate', 'ditch', 'leigh', 'sey', 'stick', 'ering', 'ington']

prefix = np.random.choice(Uk_pre)
sufix = np.random.choice(Uk_suf)

print(prefix + sufix)

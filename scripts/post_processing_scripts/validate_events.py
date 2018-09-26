art_event_file = 'events.list'
checked_event_file = 'checked_events.list'

with open (art_event_file) as f:
    event_list = []
    for line in f:
        init_min = line.split()[0] 
        saddle = line.split()[1] 
        final_min = line.split()[2] 
        event_list.append((init_min, saddle, final_min))

event_string = ''
for (init_min, saddle, final_min) in event_list:
    init_good = False
    saddle_good = False
    final_good = False 
    event_good = False
    with open (init_min + '_results.txt') as ini:
        freq_verdict = ''
        opt_verdict = ''
        for line in ini:
            if 'Frequency verdict' in line:
                freq_verdict = line.split(':')[1]
            if 'Optimization verdict' in line:
                opt_verdict = line.split(':')[1]
            if 'Indeed' in freq_verdict and 'is OK' in opt_verdict: 
                init_good = True

    with open (saddle + '_results.txt') as sad:
        freq_verdict = ''
        opt_verdict = ''
        for line in sad:
            if 'Frequency verdict' in line:
                freq_verdict = line.split(':')[1]
            if 'Optimization verdict' in line:
                opt_verdict = line.split(':')[1]
            if 'Indeed' in freq_verdict and 'is OK' in opt_verdict: 
                saddle_good = True

    with open (final_min + '_results.txt') as final:
        freq_verdict = ''
        opt_verdict = ''
        for line in final:
            if 'Frequency verdict' in line:
                freq_verdict = line.split(':')[1]
            if 'Optimization verdict' in line:
                opt_verdict = line.split(':')[1]
            if 'Indeed' in freq_verdict and 'is OK' in opt_verdict: 
                final_good = True

            if(init_good) and (saddle_good) and (final_good):
                event_good = True

    if (event_good):
        event_string += init_min + '\t' + saddle + '\t' + final_min + '\n' 

with open(checked_event_file, 'w+') as chk:
	for line in event_string.splitlines(True):
    		chk.write(line)


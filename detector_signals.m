function detected_signals = detector_signals(freq_axis, amplitude_spectrum, window_size, threshold_db)
    power_spectrum = amplitude_spectrum.^2;
    spectrum_db = 10*log10(power_spectrum + eps);
    noise_floor = median(spectrum_db);
    detection_threshold = noise_floor + threshold_db * 0.9; 
    
    detected_signals = struct('center_freq', {}, 'bandwidth', {}, 'start_freq', {}, 'end_freq', {}, ...
        'peak_power', {}, 'snr', {}, 'type', {}, 'is_nested', {}, 'parent_center', {});
                        
    min_peak_dist = max(5, round(window_size/3)); 
    [peaks, locs] = findpeaks(spectrum_db, 'MinPeakHeight', detection_threshold, 'MinPeakDistance', min_peak_dist);
    
    for i = 1:length(locs)
        sig = analyze_peak(spectrum_db, freq_axis, locs(i), peaks(i), window_size, noise_floor, threshold_db);
        if ~isempty(sig)
            detected_signals(end+1) = sig;
        end
    end
    
    bb_signals = find_broadband_signals(spectrum_db, freq_axis, window_size, noise_floor, threshold_db*0.9);
    detected_signals = [detected_signals, bb_signals];
    
    
    detected_signals = unify_narrow_and_broad(detected_signals, freq_axis);
    
    unique_new = struct('center_freq', {}, 'bandwidth', {}, 'start_freq', {}, 'end_freq', {}, ...
                       'peak_power', {}, 'snr', {}, 'type', {}, 'is_nested', {}, 'parent_center', {});
    
    while true
        new_narrow = find_nested_narrow_signals(detected_signals, spectrum_db, freq_axis, window_size, threshold_db*0.9, noise_floor);
        new_broad = find_nested_broad_signals(detected_signals, spectrum_db, freq_axis, window_size, threshold_db*0.9, noise_floor);
        new_sigs = [new_narrow, new_broad];
        
        len_before = length(unique_new);
        
        for k = 1:length(new_sigs)
            ns = new_sigs(k);
            if ~is_signal_present(detected_signals, ns) && ns.peak_power > noise_floor + threshold_db*0.7
                unique_new(end+1) = ns;
            end
        end
        
        if length(unique_new) == len_before
            break;
        end
        
        detected_signals = [detected_signals, unique_new(len_before+1:end)];
        detected_signals = merge_close_signals(detected_signals, freq_axis);
    end
    
    detected_signals = unify_narrow_and_broad(detected_signals, freq_axis);
    detected_signals = correct_broadband_centers(detected_signals, spectrum_db, freq_axis);
    detected_signals = merge_broadband_signals(detected_signals);
    detected_signals = mark_nested_signals_extended(detected_signals);
    detected_signals = remove_duplicates_vectorized(detected_signals);
    
    function sig = analyze_peak(s_db, f_axis, peak_loc, peak_pwr, win_size, noise_lvl, thresh_db)
        sig = [];
        n_spec = length(s_db);
        half_win = round(win_size/2);
        start_idx = max(1, peak_loc - half_win);
        end_idx = min(n_spec, peak_loc + half_win);
        region = s_db(start_idx:end_idx);
        freqs = f_axis(start_idx:end_idx);
        level_3db = peak_pwr - 3;
        above_3db = region >= level_3db;
        if ~any(above_3db)
            return;
        end
        start_rel = find(above_3db, 1, 'first');
        end_rel = find(above_3db, 1, 'last');
        start_freq = freqs(start_rel);
        end_freq = freqs(end_rel);
        bw = end_freq - start_freq;
        mask = (freqs >= start_freq) & (freqs <= end_freq);
        region_for_center = region(mask);
        freqs_for_center = freqs(mask);
        center_freq = sum(freqs_for_center .* 10.^(region_for_center/10)) / sum(10.^(region_for_center/10));
        local_noise = median(region);
        if local_noise < noise_lvl
            local_noise = noise_lvl;
        end
        snr = peak_pwr - local_noise;
        min_bw = 5*(f_axis(2)-f_axis(1));  
        max_narrow_bw = (max(f_axis) - min(f_axis))/4;
        min_snr = thresh_db * 1.0;
        if bw < min_bw || snr < min_snr
            return;
        end
        if bw <= max_narrow_bw
            sig_type = 'narrowband';
        else
            sig_type = 'broadband';
        end
        sig = struct('center_freq', center_freq, 'bandwidth', bw, ...
            'start_freq', start_freq, 'end_freq', end_freq, ...
            'peak_power', peak_pwr, 'snr', snr, 'type', sig_type, ...
            'is_nested', false, 'parent_center', NaN);
    end

    function bb_sigs = find_broadband_signals(s_db, f_axis, win_size, noise_lvl, thresh_db)
        bb_sigs = struct('center_freq', {}, 'bandwidth', {}, 'start_freq', {}, 'end_freq', {}, ...
            'peak_power', {}, 'snr', {}, 'type', {}, 'is_nested', {}, 'parent_center', {});
        min_bb_bw = (max(f_axis) - min(f_axis))/12;
        mask = s_db > (noise_lvl + thresh_db*0.75);
        mask = bwareaopen(mask, win_size*2);
        cc = bwconncomp(mask);
        for i = 1:cc.NumObjects
            idxs = cc.PixelIdxList{i};
            if length(idxs) < win_size*2
                continue;
            end
            start_idx = idxs(1);
            end_idx = idxs(end);
            start_freq = f_axis(start_idx);
            end_freq = f_axis(end_idx);
            bw = end_freq - start_freq;
            if bw < min_bb_bw
                continue;
            end
            region_db = s_db(start_idx:end_idx);
            peak_power = max(region_db);
            snr = peak_power - noise_lvl;
            if snr < thresh_db*0.75
                continue;
            end
            center_freq = (start_freq + end_freq)/2;
            bb_sigs(end+1) = struct('center_freq', center_freq, 'bandwidth', bw, ...
                'start_freq', start_freq, 'end_freq', end_freq, ...
                'peak_power', peak_power, 'snr', snr, ...
                'type', 'broadband', 'is_nested', false, 'parent_center', NaN);
        end
    end
    
    function signals_out = merge_close_signals(signals_in, f_axis)
        if isempty(signals_in)
            signals_out = signals_in;
            return;
        end
        [~, idx_sort] = sort([signals_in.start_freq]);
        signals_in = signals_in(idx_sort);
        to_remove = false(1,length(signals_in));
        freq_step = f_axis(2) - f_axis(1);
        freq_tol_merge = 3 * freq_step;
        for i = 1:length(signals_in)-1
            if to_remove(i)
                continue;
            end
            sig_i = signals_in(i);
            for j = i+1:length(signals_in)
                if to_remove(j)
                    continue;
                end
                sig_j = signals_in(j);
                if strcmp(sig_i.type,'narrowband') && strcmp(sig_j.type,'narrowband')
                    if (sig_j.start_freq - sig_i.end_freq) <= freq_tol_merge
                        new_start = min(sig_i.start_freq, sig_j.start_freq);
                        new_end = max(sig_i.end_freq, sig_j.end_freq);
                        new_bw = new_end - new_start;
                        new_center = (sig_i.center_freq*sig_i.bandwidth + sig_j.center_freq*sig_j.bandwidth)/(sig_i.bandwidth+sig_j.bandwidth);
                        new_peak = max(sig_i.peak_power, sig_j.peak_power);
                        new_snr = max(sig_i.snr, sig_j.snr);
                        signals_in(i).start_freq = new_start;
                        signals_in(i).end_freq = new_end;
                        signals_in(i).bandwidth = new_bw;
                        signals_in(i).center_freq = new_center;
                        signals_in(i).peak_power = new_peak;
                        signals_in(i).snr = new_snr;
                        to_remove(j) = true;
                    else
                        break;
                    end
                end
            end
        end
        signals_out = signals_in(~to_remove);
    end
    
    function signals = unify_narrow_and_broad(signals, f_axis)
        bb_mask = strcmp({signals.type}, 'broadband');
        nb_mask = strcmp({signals.type}, 'narrowband') & ~[signals.is_nested];
        bb_signals = signals(bb_mask);
        nb_signals = signals(nb_mask);
        freq_tol = (f_axis(2) - f_axis(1)) * 3.5;
        bw_tol = freq_tol * 3.5;
        to_remove_nb = false(1,length(nb_signals));
        for i = 1:length(nb_signals)
            nb = nb_signals(i);
            candidates_idx = find(arrayfun(@(bb) nb.center_freq >= bb.start_freq && nb.center_freq <= bb.end_freq, bb_signals));
            if isempty(candidates_idx)
                continue;
            end
            [~, idx_min] = min(abs([bb_signals(candidates_idx).center_freq] - nb.center_freq));
            bb_idx = candidates_idx(idx_min);
            bb = bb_signals(bb_idx);
            center_diff = abs(nb.center_freq - bb.center_freq);
            bw_diff = abs(nb.bandwidth - bb.bandwidth);
            if center_diff <= freq_tol && bw_diff <= bw_tol
                new_center = (bb.center_freq*bb.bandwidth + nb.center_freq*nb.bandwidth) / (bb.bandwidth + nb.bandwidth);
                new_bw = max(bb.bandwidth, nb.bandwidth);
                new_start = min(bb.start_freq, nb.start_freq);
                new_end = max(bb.end_freq, nb.end_freq);
                new_peak = max(bb.peak_power, nb.peak_power);
                new_snr = max(bb.snr, nb.snr);
                bb_signals(bb_idx).center_freq = new_center;
                bb_signals(bb_idx).bandwidth = new_bw;
                bb_signals(bb_idx).start_freq = new_start;
                bb_signals(bb_idx).end_freq = new_end;
                bb_signals(bb_idx).peak_power = new_peak;
                bb_signals(bb_idx).snr = new_snr;
                to_remove_nb(i) = true;
            end
        end
        nb_signals(to_remove_nb) = [];
        signals = [bb_signals, nb_signals, signals(~bb_mask & ~nb_mask)];
    end
    
    function new_signals = find_nested_narrow_signals(signals, s_db, f_axis, win_size, thresh_db, noise_level)
        bb_signals_in_list = signals(strcmp({signals.type}, 'broadband'));
        new_signals = repmat(struct('center_freq', NaN, 'bandwidth', NaN, 'start_freq', NaN, ...
            'end_freq', NaN, 'peak_power', NaN, 'snr', NaN, 'type', '', ...
            'is_nested', false, 'parent_center', NaN), 0);
        for i = 1:length(bb_signals_in_list)
            bb = bb_signals_in_list(i);
            idx_start = find(f_axis >= bb.start_freq, 1, 'first');
            idx_end = find(f_axis <= bb.end_freq, 1, 'last');
            if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
                continue;
            end
            bb_spec = s_db(idx_start:idx_end);
            bb_freqs = f_axis(idx_start:idx_end);
            med_win_for_bg = max(round(length(bb_spec)/15), 3);
            local_bg = movmedian(bb_spec, med_win_for_bg);
            norm_spec = bb_spec - local_bg;
            local_noise_nested = median(norm_spec(norm_spec < 0));
            peak_thresh_nested = max(thresh_db*0.85, local_noise_nested + thresh_db*0.75);
            min_dist_nested = max(5, round(win_size/3));
            [peaks, locs] = findpeaks(norm_spec, 'MinPeakHeight', peak_thresh_nested, 'MinPeakDistance', min_dist_nested);
            for j = 1:length(locs)
                nested_sig = analyze_peak(norm_spec, bb_freqs, locs(j), peaks(j), max(round(win_size/3),3), local_noise_nested, thresh_db*0.85);
                if ~isempty(nested_sig) ...
                    && strcmp(nested_sig.type, 'narrowband') ...
                    && nested_sig.bandwidth < 0.8*bb.bandwidth ...
                    && nested_sig.bandwidth > (f_axis(2)-f_axis(1)) ...
                    && nested_sig.snr > thresh_db*0.75 ...
                    && nested_sig.peak_power > noise_level + thresh_db*0.75
                    if nested_sig.start_freq >= bb.start_freq && nested_sig.end_freq <= bb.end_freq
                        if ~is_signal_present(new_signals, nested_sig)
                            nested_sig.is_nested = true;
                            nested_sig.parent_center = bb.center_freq;
                            new_signals(end+1) = nested_sig;
                        end
                    end
                end
            end
        end
    end
    
    function new_signals = find_nested_broad_signals(signals, s_db, f_axis, win_size, thresh_db, noise_level)
        bb_signals_in_list = signals(strcmp({signals.type}, 'broadband'));
        new_signals = repmat(struct('center_freq', NaN, 'bandwidth', NaN, 'start_freq', NaN, ...
            'end_freq', NaN, 'peak_power', NaN, 'snr', NaN, 'type', '', ...
            'is_nested', false, 'parent_center', NaN), 0);
        min_nested_bw_ratio = 0.3;
        max_nested_bw_ratio = 0.99; 
        
        max_width_diff_hz = 5; 
        
        for i = 1:length(bb_signals_in_list)
            bb = bb_signals_in_list(i);
            idx_start = find(f_axis >= bb.start_freq, 1, 'first');
            idx_end = find(f_axis <= bb.end_freq, 1, 'last');
            if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
                continue;
            end
            spec = s_db(idx_start:idx_end);
            freqs = f_axis(idx_start:idx_end);
            noise_est = median(spec);
            threshold = noise_est + thresh_db*0.75;
            mask = spec > threshold;
            mask = bwareaopen(mask, round(win_size*2));
            cc = bwconncomp(mask);
            for cc_i = 1:cc.NumObjects
                pix = cc.PixelIdxList{cc_i};
                sub_bw = freqs(pix(end)) - freqs(pix(1));
                width_diff = abs(sub_bw - bb.bandwidth);
                if sub_bw >= min_nested_bw_ratio * bb.bandwidth && ...
                   (sub_bw <= max_nested_bw_ratio * bb.bandwidth || width_diff <= max_width_diff_hz)
                    sub_start = freqs(pix(1));
                    sub_end = freqs(pix(end));
                    sub_center = (sub_start + sub_end)/2;
                    sub_power = max(spec(pix));
                    sub_snr = sub_power - noise_level;
                    if sub_snr >= thresh_db*0.5 && sub_start >= bb.start_freq && sub_end <= bb.end_freq
                        new_sig = struct('center_freq', sub_center, 'bandwidth', sub_bw, ...
                            'start_freq', sub_start, 'end_freq', sub_end, ...
                            'peak_power', sub_power, 'snr', sub_snr, ...
                            'type', 'broadband', 'is_nested', true, 'parent_center', bb.center_freq);
                        if ~is_signal_present(new_signals, new_sig)
                            new_signals(end+1) = new_sig;
                        end
                    end
                end
            end
        end
    end
    
    function present = is_signal_present(signal_array, signal)
        present = false;
        if isempty(signal_array)
            return;
        end
        freq_tol_check = (freq_axis(2) - freq_axis(1)) * 0.05;
        bw_tol_check = freq_tol_check;
        for k = 1:length(signal_array)
            if abs(signal_array(k).center_freq - signal.center_freq) <= freq_tol_check && ...
               abs(signal_array(k).bandwidth - signal.bandwidth) <= bw_tol_check
                present = true;
                return;
            end
        end
    end
    
    function signals = correct_broadband_centers(signals, s_db, f_axis)
        for i = 1:length(signals)
            if ~strcmp(signals(i).type, 'broadband')
                continue;
            end
            idx_start = find(f_axis >= signals(i).start_freq, 1, 'first');
            idx_end = find(f_axis <= signals(i).end_freq, 1, 'last');
            if isempty(idx_start) || isempty(idx_end) || idx_start >= idx_end
                continue;
            end
            broad_spec = s_db(idx_start:idx_end);
            broad_freqs = f_axis(idx_start:idx_end);
            nested_sigs = signals([signals.is_nested] & ...
                ([signals.start_freq] >= signals(i).start_freq) & ...
                ([signals.end_freq] <= signals(i).end_freq));
            cleaned_spec = broad_spec;
            for k = 1:length(nested_sigs)
                ns = nested_sigs(k);
                ns_start_idx = find(broad_freqs >= ns.start_freq, 1, 'first');
                ns_end_idx = find(broad_freqs <= ns.end_freq, 1, 'last');
                if isempty(ns_start_idx) || isempty(ns_end_idx)
                    continue;
                end
                left_val = cleaned_spec(max(ns_start_idx-1,1));
                right_val = cleaned_spec(min(ns_end_idx+1,length(cleaned_spec)));
                interp_vals = linspace(left_val, right_val, ns_end_idx - ns_start_idx + 1);
                cleaned_spec(ns_start_idx:ns_end_idx) = interp_vals;
            end
            signals(i).center_freq = (signals(i).start_freq + signals(i).end_freq) / 2;
        end
    end
    
    function signals = merge_broadband_signals(signals)
        nested_mask = [signals.is_nested];
        non_nested = signals(~nested_mask);
        nested = signals(nested_mask);
        bb_mask = strcmp({non_nested.type},'broadband');
        bb = non_nested(bb_mask);
        nb = non_nested(~bb_mask);
        if isempty(bb)
            merged_bb = [];
        else
            [~, idx] = sort([bb.start_freq]);
            bb = bb(idx);
            merged_bb = bb(1);
            for i = 2:length(bb)
                last = merged_bb(end);
                curr = bb(i);
                if curr.start_freq <= last.end_freq
                    new_start = min(last.start_freq, curr.start_freq);
                    new_end = max(last.end_freq, curr.end_freq);
                    new_center = (last.center_freq*last.bandwidth + curr.center_freq*curr.bandwidth) / (last.bandwidth + curr.bandwidth);
                    new_pwr = max(last.peak_power, curr.peak_power);
                    new_snr = max(last.snr, curr.snr);
                    merged_bb(end) = struct('center_freq', new_center, 'bandwidth', new_end-new_start, ...
                        'start_freq', new_start, 'end_freq', new_end, ...
                        'peak_power', new_pwr, 'snr', new_snr, 'type', 'broadband', ...
                        'is_nested', false, 'parent_center', NaN);
                else
                    merged_bb(end+1) = curr;
                end
            end
        end
        signals = [merged_bb, nb, nested];
    end
    
    function signals = mark_nested_signals_extended(signals)
        n = length(signals);
        if n < 2
            return;
        end
        [~, idx_sort] = sort([signals.center_freq]);
        signals = signals(idx_sort);
        for i = 1:n-1
            if signals(i).is_nested
                continue;
            end
            bw_i = signals(i).bandwidth;
            start_i = signals(i).start_freq;
            end_i = signals(i).end_freq;
            for j = i+1:n
                if signals(j).is_nested
                    continue;
                end
                if signals(j).start_freq > end_i
                    break;
                end
                if (signals(j).start_freq >= start_i + 0.03*bw_i) && ...
                   (signals(j).end_freq <= end_i - 0.03*bw_i) && ...
                   (signals(j).bandwidth < 0.95*bw_i)
                    signals(j).is_nested = true;
                    signals(j).parent_center = signals(i).center_freq;
                end
            end
        end
    end
    
    function signals_out = remove_duplicates_vectorized(signals_in)
        if isempty(signals_in)
            signals_out = signals_in;
            return;
        end
        center_freqs = [signals_in.center_freq];
        bandwidths = [signals_in.bandwidth];
        types = {signals_in.type};
        freq_tol = 0.05;
        bw_tol = 0.05;
        keep_mask = true(1, numel(signals_in));
        for i = 1:numel(signals_in)
            if ~keep_mask(i)
                continue;
            end
            same_type_idx = find(strcmp(types, types{i}));
            same_type_idx = same_type_idx(same_type_idx > i);
            for j = same_type_idx
                if ~keep_mask(j)
                    continue;
                end
                if abs(center_freqs(i) - center_freqs(j)) <= freq_tol && ...
                   abs(bandwidths(i) - bandwidths(j)) <= bw_tol
                    keep_mask(j) = false;
                end
            end
        end
        signals_out = signals_in(keep_mask);
    end
end

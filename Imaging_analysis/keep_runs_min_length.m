function cleaned = keep_runs_min_length(mask, minlen)
    mask = mask(:)'; % force row vector
    z = [0 mask 0];      % pad with zeros at both ends
    d = diff(z);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;
    run_lengths = run_ends - run_starts + 1;
    % Keep runs at least minlen long
    ok = run_lengths >= minlen;
    cleaned = false(size(mask));
    for i = find(ok)
        cleaned(run_starts(i):run_ends(i)) = true;
    end
    cleaned = cleaned(:); % column vector output
end
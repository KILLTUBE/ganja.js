equals(simplify('e0', 0, 2, 1), 'e0');
equals(simplify('e1', 0, 2, 1), 'e1');
equals(simplify('e2', 0, 2, 1), 'e2');
equals(simplify('e01', 0, 2, 1), 'e01');
equals(simplify('e02', 0, 2, 1), 'e02');
equals(simplify('e12', 0, 2, 1), 'e12');
equals(simplify('e012', 0, 2, 1), 'e012');
equals(simplify('e12e0', undefined, undefined, undefined), 'e012');
equals(simplify('e02e1', undefined, undefined, undefined), '-e012');
equals(simplify('e01e2', undefined, undefined, undefined), 'e012');
equals(simplify('e2e01', undefined, undefined, undefined), 'e012');
equals(simplify('e1e02', undefined, undefined, undefined), '-e012');
equals(simplify('e0e12', undefined, undefined, undefined), 'e012');
equals(simplify('11', 0, 2, 1), '-1');
equals(simplify('e0e0', 0, 2, 1), '0');
equals(simplify('e1e1', 0, 2, 1), '-1');
equals(simplify('e2e2', 0, 2, 1), '-1');
equals(simplify('e01e01', 0, 2, 1), '0');
equals(simplify('e02e02', 0, 2, 1), '0');
equals(simplify('e12e12', 0, 2, 1), '-1');
equals(simplify('e012e012', 0, 2, 1), '0');
equals(simplify('e0e0', 0, 2, 1), '0');
equals(simplify('e0e1', 0, 2, 1), 'e01');
equals(simplify('e0e2', 0, 2, 1), 'e02');
equals(simplify('e0e01', 0, 2, 1), '0');
equals(simplify('e0e02', 0, 2, 1), '0');
equals(simplify('e0e12', 0, 2, 1), 'e012');
equals(simplify('e0e012', 0, 2, 1), '0');
equals(simplify('e1e0', 0, 2, 1), '-e01');
equals(simplify('e1e1', 0, 2, 1), '-1');
equals(simplify('e1e2', 0, 2, 1), 'e12');
equals(simplify('e1e01', 0, 2, 1), 'e0');
equals(simplify('e1e02', 0, 2, 1), '-e012');
equals(simplify('e1e12', 0, 2, 1), '-e2');
equals(simplify('e1e012', 0, 2, 1), 'e02');
equals(simplify('e2e0', 0, 2, 1), '-e02');
equals(simplify('e2e1', 0, 2, 1), '-e12');
equals(simplify('e2e2', 0, 2, 1), '-1');
equals(simplify('e2e01', 0, 2, 1), 'e012');
equals(simplify('e2e02', 0, 2, 1), 'e0');
equals(simplify('e2e12', 0, 2, 1), 'e1');
equals(simplify('e2e012', 0, 2, 1), '-e01');
equals(simplify('e01e0', 0, 2, 1), '0');
equals(simplify('e01e1', 0, 2, 1), '-e0');
equals(simplify('e01e2', 0, 2, 1), 'e012');
equals(simplify('e01e01', 0, 2, 1), '0');
equals(simplify('e01e02', 0, 2, 1), '0');
equals(simplify('e01e12', 0, 2, 1), '-e02');
equals(simplify('e01e012', 0, 2, 1), '0');
equals(simplify('e02e0', 0, 2, 1), '0');
equals(simplify('e02e1', 0, 2, 1), '-e012');
equals(simplify('e02e2', 0, 2, 1), '-e0');
equals(simplify('e02e01', 0, 2, 1), '0');
equals(simplify('e02e02', 0, 2, 1), '0');
equals(simplify('e02e12', 0, 2, 1), 'e01');
equals(simplify('e02e012', 0, 2, 1), '0');
equals(simplify('e12e0', 0, 2, 1), 'e012');
equals(simplify('e12e1', 0, 2, 1), 'e2');
equals(simplify('e12e2', 0, 2, 1), '-e1');
equals(simplify('e12e01', 0, 2, 1), 'e02');
equals(simplify('e12e02', 0, 2, 1), '-e01');
equals(simplify('e12e12', 0, 2, 1), '-1');
equals(simplify('e12e012', 0, 2, 1), '-e0');
equals(simplify('e012e0', 0, 2, 1), '0');
equals(simplify('e012e1', 0, 2, 1), 'e02');
equals(simplify('e012e2', 0, 2, 1), '-e01');
equals(simplify('e012e01', 0, 2, 1), '0');
equals(simplify('e012e02', 0, 2, 1), '0');
equals(simplify('e012e12', 0, 2, 1), '-e0');
equals(simplify('e012e012', 0, 2, 1), '0');

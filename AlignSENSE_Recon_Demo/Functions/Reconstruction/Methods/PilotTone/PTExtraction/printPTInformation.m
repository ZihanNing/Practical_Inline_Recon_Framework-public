

function [] = printPTInformation(rec)

fprintf('Pilot Tone signal detected:\n');
fprintf('      idxRO:%.0f\n', rec.PT.idxRO);
fprintf('      factorFOV:%.2f\n', rec.PT.factorFOV);
fprintf('      FWHMHz:%.2f\n', rec.PT.FWHMHz);
fprintf('      isAveragedMB:%.0f\n', rec.PT.isAveragedMB);
fprintf('      isRelativePhase:%.0f\n', rec.PT.isRelativePhase);
fprintf('      relativePhaseCoilIdx:%.0f\n', rec.PT.relativePhaseCoilIdx);
fprintf('\n\n');

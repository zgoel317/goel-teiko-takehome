-- Normalized, long-format schema for immune cell population data.
-- Foreign keys are enabled per-connection in db/connection.py (a PRAGMA in a
-- schema file does not persist), so this DDL relies on that being set.

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT,
    treatment  TEXT,
    response   TEXT              -- nullable: healthy / 'none' treatment have no response
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER
);

CREATE TABLE populations (
    population_id INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id     TEXT    NOT NULL REFERENCES samples(sample_id),
    population_id INTEGER NOT NULL REFERENCES populations(population_id),
    count         INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_project   ON subjects(project_id);
CREATE INDEX idx_subjects_condition ON subjects(condition);
CREATE INDEX idx_subjects_treatment ON subjects(treatment);
CREATE INDEX idx_subjects_response  ON subjects(response);
CREATE INDEX idx_samples_subject    ON samples(subject_id);
CREATE INDEX idx_samples_type       ON samples(sample_type);
CREATE INDEX idx_samples_time       ON samples(time_from_treatment_start);
CREATE INDEX idx_cell_counts_pop    ON cell_counts(population_id, sample_id);

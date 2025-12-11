"""
MSNoise Web Administration Interface - Consolidated Version

This module provides a Flask-Admin based web interface for managing MSNoise
configurations, workflow chains, filters, and data processing jobs.

Key Features:
- Unified Configuration Management (global and workflow-specific)
- Dynamic Workflow Chain Creation and Management
- Enhanced Filter Management with Workflow Associations
- Simplified Base Classes with Modern Flask-Admin
- Generic Workflow Handling driven by Database Schema
"""

import datetime
import logging
import os
import sys
import traceback
from collections import OrderedDict

from flask import Flask, request, redirect, url_for, flash, jsonify, render_template_string
from flask_admin import Admin, BaseView, expose, AdminIndexView
from flask_admin.contrib.sqla import ModelView
from flask_admin.actions import action
from flask_admin.form import Select2Widget
from flask_admin.model import typefmt
from flask_admin import helpers as admin_helpers
from flask_babel import Babel, gettext, ngettext, lazy_gettext
from markupsafe import Markup
SafeMarkup = Markup
from sqlalchemy import and_, or_, func, desc, asc
from sqlalchemy.orm import sessionmaker
from wtforms import Form, StringField, TextAreaField, SelectField, BooleanField, IntegerField, FloatField
from wtforms.validators import DataRequired, Optional, NumberRange, ValidationError
from wtforms.widgets import TextArea, Select, CheckboxInput, NumberInput

from .api import connect, get_config, get_logger
from .api import get_config_categories_definition
from .msnoise_table_def import declare_tables


# Initialize Flask app and extensions
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'msnoise-admin-secret-key')
app.config['FLASK_ADMIN_SWATCH'] = 'cerulean'

# Babel for internationalization
babel = Babel(app)

# Database connection
db_session = connect()

# Get schema from declare_tables
schema = declare_tables()

# Logger
logger = get_logger('msnoise.admin', loglevel="INFO")

@app.route('/admin/job-stats')
def job_stats():
    """API endpoint for job statistics"""
    try:
        from .api import get_workflow_job_counts
        db = connect()

        # Get all workflows
        workflows = db.query(Job.workflow_id).distinct().all()
        stats = {}

        for workflow in workflows:
            workflow_id = workflow[0]
            stats[workflow_id] = get_workflow_job_counts(db, workflow_id)

        return jsonify(stats)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/admin/job-dependencies/<int:job_id>')
def job_dependencies(job_id):
    """API endpoint for job dependencies"""
    try:
        db = connect()
        job = db.query(Job).filter(Job.ref == job_id).first()

        if not job:
            return jsonify({'error': 'Job not found'}), 404

        dependencies = {
            'predecessors': [],
            'successors': []
        }

        if job.workflow_step:
            # Get predecessor jobs
            for link in job.workflow_step.incoming_links:
                if link.is_active:
                    pred_jobs = db.query(Job).filter(
                        Job.day == job.day,
                        Job.pair == job.pair,
                        Job.step_id == link.from_step_id,
                        Job.workflow_id == job.workflow_id
                    ).all()

                    for pred_job in pred_jobs:
                        dependencies['predecessors'].append({
                            'id': pred_job.ref,
                            'step_name': pred_job.workflow_step.step_name if pred_job.workflow_step else 'Unknown',
                            'flag': pred_job.flag
                        })

            # Get successor jobs
            for link in job.workflow_step.outgoing_links:
                if link.is_active:
                    succ_jobs = db.query(Job).filter(
                        Job.day == job.day,
                        Job.pair == job.pair,
                        Job.step_id == link.to_step_id,
                        Job.workflow_id == job.workflow_id
                    ).all()

                    for succ_job in succ_jobs:
                        dependencies['successors'].append({
                            'id': succ_job.ref,
                            'step_name': succ_job.workflow_step.step_name if succ_job.workflow_step else 'Unknown',
                            'flag': succ_job.flag
                        })

        return jsonify(dependencies)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# ============================================================================
# Base Classes and Utilities
# ============================================================================

class BaseModelView(ModelView):
    """
    Enhanced base model view with common functionality for all MSNoise models.
    Simplified from the original complex implementation.
    """

    # Common settings
    can_export = True
    can_view_details = True
    page_size = 50
    can_set_page_size = True

    # Default column formatters
    column_formatters = {
        'used': lambda view, context, model, name: SafeMarkup('✓' if model.used else '✗') if hasattr(model, 'used') else '',
    }

    def __init__(self, model, session, *args, **kwargs):
        super(BaseModelView, self).__init__(model, session, *args, **kwargs)
        self.session = session

    def get_query(self):
        """Override to add common filters"""
        return self.session.query(self.model)

    def get_count_query(self):
        """Override to add common filters"""
        return self.session.query(func.count('*')).select_from(self.model)


class MSNoiseAdminIndexView(AdminIndexView):
    """Custom admin index view with MSNoise-specific dashboard"""

    @expose('/')
    def index(self):
        """Display MSNoise dashboard"""
        # Get basic statistics
        stats = {}
        try:
            stats['stations'] = db_session.query(schema.Station).filter(schema.Station.used == True).count()
            # stats['filters'] = db_session.query(schema.Filter).filter(schema.Filter.used == True).count()
            # stats['workflow_steps'] = db_session.query(schema.WorkflowSteps).filter(schema.WorkflowSteps.used == True).count()
            stats['jobs_total'] = db_session.query(schema.Job).count()
            stats['jobs_todo'] = db_session.query(schema.Job).filter(schema.Job.flag == 'T').count()
            stats['jobs_done'] = db_session.query(schema.Job).filter(schema.Job.flag == 'D').count()
            stats['config_params'] = db_session.query(schema.Config).count()
        except Exception as e:
            logger.error(f"Error getting dashboard stats: {e}")
            stats = {}

        return self.render('admin/index.html', stats=stats)


# ============================================================================
# 1. Single ConfigView - Unified Configuration Management
# ============================================================================

class ConfigView(BaseModelView):
    """
    Unified configuration management for both global and workflow-specific settings.
    This replaces all the separate DvvMwcsView, DvvStretchingView, etc. classes.
    """

    column_list = ['category', 'name', 'set_number', 'value', 'description']
    column_searchable_list = ['name', 'category', 'value', 'description']
    column_filters = ['category', 'set_number']
    column_sortable_list = ['name', 'category', 'set_number']

    form_columns = ['name', 'category', 'set_number', 'value', 'param_type',
                    'default_value', 'description', 'units', 'possible_values']

    column_descriptions = {
        'name': 'Parameter name (e.g., maxlag, dtt_minlag)',
        'category': 'Configuration category (global, mwcs, stretching, etc.)',
        'set_number': 'Configuration set number (NULL for global configs)',
        'value': 'Current parameter value',
        'param_type': 'Data type of the parameter',
        'description': 'Description of what this parameter does',
    }

    form_choices = {
        'category': get_config_categories_definition(),
        'param_type': [
            ('str', 'String'),
            ('int', 'Integer'),
            ('float', 'Float'),
            ('bool', 'Boolean'),
        ]
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'value': lambda view, context, model, name: SafeMarkup(f'<code>{model.value}</code>') if model.value else '',
        'category': lambda view, context, model, name: SafeMarkup(f'<span class="label label-info">{model.category}</span>'),
        'param_type': lambda view, context, model, name: SafeMarkup(f'<span class="label label-default">{model.param_type}</span>'),
        'set_number': lambda view, context, model, name: str(model.set_number) if model.set_number is not None else 'Global',
    })

    def __init__(self, session, *args, **kwargs):
        super(ConfigView, self).__init__(schema.Config, session, *args, **kwargs)

    def create_model(self, form):
        """Create new configuration parameter with validation"""
        try:
            model = self.model()
            form.populate_obj(model)

            # Validate value type
            if model.param_type == 'int':
                int(model.value)
            elif model.param_type == 'float':
                float(model.value)
            elif model.param_type == 'bool':
                model.value = str(model.value).lower() in ('true', 'yes', '1', 'on')

            # Validate against possible values if defined
            model.validate_value()

            self.session.add(model)
            self.session.commit()
            return model
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to create configuration: {str(e)}', 'error')
            return False

    def update_model(self, form, model):
        """Update configuration parameter with validation"""
        try:
            form.populate_obj(model)

            # Validate value type
            if model.param_type == 'int':
                int(model.value)
            elif model.param_type == 'float':
                float(model.value)
            elif model.param_type == 'bool':
                model.value = str(model.value).lower() in ('true', 'yes', '1', 'on')

            # Validate against possible values if defined
            model.validate_value()

            self.session.commit()
            return True
        except Exception as e:
            self.session.rollback()
            flash(f'Failed to update configuration: {str(e)}', 'error')
            return False

    @expose('/bulk_update/', methods=['GET', 'POST'])
    def bulk_update(self):
        """Bulk update configuration parameters"""
        if request.method == 'POST':
            try:
                updates = request.form.to_dict()
                for key, value in updates.items():
                    if key.startswith('config_'):
                        config_id = int(key.split('_')[1])
                        config = self.session.query(schema.Config).get(config_id)
                        if config:
                            config.value = value
                            config.updated_at = datetime.datetime.utcnow()

                self.session.commit()
                flash('Configuration updated successfully', 'success')
            except Exception as e:
                self.session.rollback()
                flash(f'Failed to update configuration: {str(e)}', 'error')

            return redirect(url_for('.index_view'))

        # GET request - show bulk update form
        configs = self.session.query(schema.Config).all()
        return self.render('admin/config_bulk_update.html', configs=configs)


# ============================================================================
# 4. Simplified Base Classes - Standard Model Views
# ============================================================================

class StationView(BaseModelView):
    """Simplified station management"""

    column_list = ['ref', 'net', 'sta', 'X', 'Y', 'altitude', 'coordinates', 'used']
    column_searchable_list = ['net', 'sta']
    column_filters = ['net', 'coordinates', 'used']
    column_sortable_list = ['ref', 'net', 'sta', 'X', 'Y', 'used']

    form_columns = ['net', 'sta', 'X', 'Y', 'altitude', 'coordinates',
                    'used_location_codes', 'used_channel_names', 'used']

    def __init__(self, session, *args, **kwargs):
        super(StationView, self).__init__(schema.Station, session, *args, **kwargs)


class JobView(BaseModelView):
    """Enhanced Job view with workflow support"""

    def __init__(self, session, *args, **kwargs):
        # Import Job model here to avoid circular imports
        from .msnoise_table_def import Job
        # Call parent with model and session
        super(JobView, self).__init__(Job, session, *args, **kwargs)

    # Update column list to include new workflow fields
    column_list = ['ref', 'day', 'pair', 'workflow_id', 'step_name', 'priority', 'jobtype', 'flag', 'lastmod']

    # Add searchable columns
    column_searchable_list = ['day', 'pair', 'workflow_id', 'jobtype']

    # Add filters for better navigation
    column_filters = ['flag', 'workflow_id', 'jobtype', 'day', 'priority']

    # Sort by priority and last modified by default
    column_default_sort = [('priority', True), ('lastmod', False)]

    # Add column labels for better readability
    column_labels = {
        'ref': 'Job ID',
        'day': 'Date',
        'pair': 'Station Pair',
        'workflow_id': 'Workflow',
        'step_name': 'Step',
        'priority': 'Priority',
        'jobtype': 'Job Type',
        'flag': 'Status',
        'lastmod': 'Last Modified'
    }

    # Add descriptions for columns
    column_descriptions = {
        'workflow_id': 'Workflow identifier this job belongs to',
        'step_name': 'Workflow step name',
        'priority': 'Job priority (higher = more important)',
        'flag': 'T=Todo, I=In Progress, D=Done'
    }

    # Custom column formatting
    column_formatters = {
        'flag': lambda v, c, m, p: {
            'T': Markup('<span class="label label-warning">Todo</span>'),
            'I': Markup('<span class="label label-info">In Progress</span>'),
            'D': Markup('<span class="label label-success">Done</span>')
        }.get(m.flag, m.flag),
        'priority': lambda v, c, m, p: Markup(f'<span class="badge badge-{"danger" if m.priority > 5 else "warning" if m.priority > 0 else "default"}">{m.priority}</span>'),
        'step_name': lambda v, c, m, p: Markup(m.workflow_step.step_name if hasattr(m, 'workflow_step') and m.workflow_step else 'N/A'),

    }

    # Form configuration
    form_columns = ['day', 'pair', 'workflow_id', 'step_id', 'priority', 'jobtype', 'flag']

    # Custom query to join with WorkflowStep
    def get_query(self):
        from .msnoise_table_def import WorkflowStep
        return self.session.query(self.model).outerjoin(WorkflowStep)

    def get_count_query(self):
        from .msnoise_table_def import WorkflowStep
        from sqlalchemy import func
        return self.session.query(func.count('*')).select_from(self.model).outerjoin(WorkflowStep)

    # Add custom actions with proper decorator
    @action('reset_jobs', 'Reset Selected Jobs', 'Are you sure you want to reset the selected jobs?')
    def action_reset_jobs(self, ids):
        try:
            from datetime import datetime
            query = self.session.query(self.model).filter(self.model.ref.in_(ids))
            count = 0
            for job in query.all():
                job.flag = 'T'
                job.lastmod = datetime.utcnow()
                count += 1
            self.session.commit()
            flash(ngettext('Job was successfully reset.',
                           '%(count)s jobs were successfully reset.',
                           count, count=count), 'success')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to reset jobs. %(error)s', error=str(ex)), 'error')

    @action('mark_done', 'Mark as Done', 'Are you sure you want to mark the selected jobs as done?')
    def action_mark_done(self, ids):
        try:
            from datetime import datetime
            query = self.session.query(self.model).filter(self.model.ref.in_(ids))
            count = 0
            for job in query.all():
                job.flag = 'D'
                job.lastmod = datetime.utcnow()
                count += 1
            self.session.commit()
            flash(ngettext('Job was successfully marked as done.',
                           '%(count)s jobs were successfully marked as done.',
                           count, count=count), 'success')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to mark jobs as done. %(error)s', error=str(ex)), 'error')

    @action('mark_todo', 'Mark as Todo', 'Are you sure you want to mark the selected jobs as todo?')
    def action_mark_todo(self, ids):
        try:
            from datetime import datetime
            query = self.session.query(self.model).filter(self.model.ref.in_(ids))
            count = 0
            for job in query.all():
                job.flag = 'T'
                job.lastmod = datetime.utcnow()
                count += 1
            self.session.commit()
            flash(ngettext('Job was successfully marked as todo.',
                           '%(count)s jobs were successfully marked as todo.',
                           count, count=count), 'success')
        except Exception as ex:
            if not self.handle_view_exception(ex):
                raise
            flash(gettext('Failed to mark jobs as todo. %(error)s', error=str(ex)), 'error')

    # Add workflow statistics to context
    def get_list_stats(self):
        """Get job statistics for the current workflow"""
        from sqlalchemy import func
        stats = {}
        workflows = self.session.query(self.model.workflow_id).distinct().all()

        for workflow in workflows:
            workflow_id = workflow[0]
            counts = self.session.query(
                self.model.flag,
                func.count(self.model.ref).label('count')
            ).filter(
                self.model.workflow_id == workflow_id
            ).group_by(self.model.flag).all()

            stats[workflow_id] = {
                'T': 0, 'I': 0, 'D': 0
            }
            for flag, count in counts:
                stats[workflow_id][flag] = count

        return stats

class DataAvailabilityView(BaseModelView):
    """Simplified data availability management"""

    column_list = ['ref', 'net', 'sta', 'loc', 'chan', 'starttime', 'endtime', 'samplerate', 'flag']
    column_searchable_list = ['net', 'sta', 'loc', 'chan', 'file']
    column_filters = ['net', 'sta', 'chan', 'flag']
    column_sortable_list = ['ref', 'net', 'sta', 'starttime', 'endtime', 'samplerate']

    # Read-only view
    can_create = False
    can_edit = False
    can_delete = False

    def __init__(self, session, *args, **kwargs):
        super(DataAvailabilityView, self).__init__(schema.DataAvailability, session, *args, **kwargs)

from flask_admin import BaseView, expose
from flask_admin.form import BaseForm
from flask import request, redirect, url_for, flash
from wtforms import Form, StringField, FloatField, IntegerField, BooleanField, SelectField, TextAreaField
from wtforms.validators import DataRequired, Optional
import datetime

class ConfigSetView(BaseView):
    """
    View for managing configuration sets as single forms.
    Allows editing all parameters in a set together.
    """

    def __init__(self, session, *args, **kwargs):
        super(ConfigSetView, self).__init__(*args, **kwargs)
        self.session = session

    @expose('/')
    def index(self):
        """Display all configuration sets grouped by category"""

        """Display all configuration sets grouped by category"""
        from .api import get_config_sets_organized

        ordered_categories = get_config_sets_organized(self.session)

        return self.render('admin/config_sets_index.html',
                           categories=ordered_categories)


    @expose('/create/<category>/', methods=['GET', 'POST'])
    def create_set(self, category):
        """Create a new configuration set for a category"""

        from .api import create_config_set
        set_number = create_config_set(self.session, category)


        if set_number is not None:
            self.session.commit()
            flash(f'Configuration set "{category}" #{set_number} created successfully', 'success')
            return redirect(url_for('.edit_set', category=category, set_number=set_number))

        return

    @expose('/view/<category>/<int:set_number>/')
    def view_set(self, category, set_number):
        """View configuration set in read-only mode"""

        # Get all configs for this category and set
        configs = self.session.query(schema.Config).filter(
            schema.Config.category == category,
            schema.Config.set_number == set_number
        ).order_by(schema.Config.name).all()

        if not configs:
            flash(f'Configuration set "{category}" #{set_number} not found', 'error')
            return redirect(url_for('.index'))

        # Group configs by logical sections if needed
        config_groups = {}
        for config in configs:
            config.description = Markup(config.description)
            # You can add logic here to group configs by prefix or type
            group_name = "Configuration Parameters"
            if group_name not in config_groups:
                config_groups[group_name] = []

            config_groups[group_name].append(config)

        return self.render('admin/config_set_view.html',
                           category=category,
                           set_number=set_number,
                           configs=configs,
                           config_groups=config_groups)


    @expose('/edit/<category>/<int:set_number>/', methods=['GET', 'POST'])
    def edit_set(self, category, set_number):
        """Edit a complete configuration set as a single form"""

        # Get all config parameters for this set
        configs = self.session.query(schema.Config).filter(
            schema.Config.category == category,
            schema.Config.set_number == set_number
        ).order_by(schema.Config.name).all()

        if not configs:
            flash(f'Configuration set "{category}" #{set_number} not found', 'error')
            return redirect(url_for('.index'))

        # Create dynamic form class
        form_class = self._create_dynamic_form(configs)
        # Initialize form with request data for POST requests
        if request.method == 'POST':
            form = form_class(request.form)
        else:
            form = form_class()

        if request.method == 'POST' and form.validate():
            try:
                # Update all config values
                for config in configs:
                    field_name = f"param_{config.ref}"
                    if hasattr(form, field_name):
                        field = getattr(form, field_name)
                        new_value = field.data

                        # Handle different field types and preserve existing values if new value is empty
                        if config.param_type == 'bool':
                            # For boolean fields
                            if isinstance(field, BooleanField):
                                # WTForms BooleanField returns True/False
                                new_value = 'Y' if new_value else 'N'
                            elif isinstance(field, SelectField):
                                # SelectField for boolean with Y/N choices
                                new_value = str(new_value) if new_value else config.value
                            else:
                                new_value = str(new_value) if new_value else config.value
                        elif config.param_type == 'int':
                            if new_value is not None and str(new_value).strip():
                                try:
                                    new_value = str(int(new_value))
                                except (ValueError, TypeError):
                                    # Keep original value if conversion fails
                                    new_value = config.value
                            else:
                                # Keep original value if empty
                                new_value = config.value
                        elif config.param_type == 'float':
                            if new_value is not None and str(new_value).strip():
                                try:
                                    new_value = str(float(new_value))
                                except (ValueError, TypeError):
                                    # Keep original value if conversion fails
                                    new_value = config.value
                            else:
                                # Keep original value if empty
                                new_value = config.value
                        else:  # string type
                            if new_value is not None:
                                new_value = str(new_value).strip()
                                # Allow empty strings for string fields
                            else:
                                new_value = config.value

                        # Only update if we have a valid new value
                        if new_value is not None:
                            config.value = new_value
                            config.updated_at = datetime.datetime.utcnow()

                            # Validate against possible values (only for non-empty values)
                            if new_value.strip():
                                try:
                                    config.validate_value()
                                except ValueError as ve:
                                    flash(f'Validation error for {config.name}: {str(ve)}', 'error')
                                    raise ve
                            # flash(f'Updated parameter "{config.name}"="{new_value}"', 'success')
                            self.session.commit()
                self.session.commit()
                flash(f'Configuration set "{category}" #{set_number} updated successfully', 'success')
                return redirect(url_for('.view_set', category=category, set_number=set_number))

            except Exception as e:
                self.session.rollback()
                flash(f'Error updating configuration: {str(e)}', 'error')

        else:
            # Populate form with current values (GET request or validation failed)
            for config in configs:
                config.description = Markup(config.description)
                field_name = f"param_{config.ref}"
                if hasattr(form, field_name):
                    field = getattr(form, field_name)
                    current_value = config.value if config.value is not None else ''
                    field.description = config.description
                    if config.param_type == 'bool':
                        # Handle boolean conversion for different formats
                        if isinstance(field, BooleanField):
                            field.data = current_value.upper() in ('Y', 'TRUE', 'YES', '1', 'ON')
                        else:  # SelectField for boolean
                            field.data = current_value
                    elif config.param_type in ['int', 'float']:
                        try:
                            if current_value and str(current_value).strip():
                                if config.param_type == 'float':
                                    field.data = float(current_value)
                                else:
                                    field.data = int(current_value)
                            else:
                                field.data = None
                        except (ValueError, TypeError):
                            # If conversion fails, set as string for display
                            field.data = current_value
                    else:
                        field.data = current_value

            # Show validation errors if any
            if request.method == 'POST':
                for field_name, errors in form.errors.items():
                    for error in errors:
                        flash(f'Validation error in {field_name}: {error}', 'error')

        return self.render('admin/config_set_edit.html',
                           form=form,
                           configs=configs,
                           category=category,
                           set_number=set_number)

    @expose('/copy/<category>/<int:set_number>/', methods=['GET', 'POST'])
    def copy_set(self, category, set_number):
        """Copy a configuration set to create a new one"""

        if request.method == 'POST':
            new_category = request.form.get('new_category', category)
            new_set_number = int(request.form.get('new_set_number', set_number + 1))

            try:
                # Check if target already exists
                existing = self.session.query(schema.Config).filter(
                    schema.Config.category == new_category,
                    schema.Config.set_number == new_set_number
                ).first()

                if existing:
                    flash(f'Configuration set "{new_category}" #{new_set_number} already exists', 'error')
                    return redirect(url_for('.index'))

                # Get source configs
                source_configs = self.session.query(schema.Config).filter(
                    schema.Config.category == category,
                    schema.Config.set_number == set_number
                ).all()

                if not source_configs:
                    flash(f'Source configuration set not found', 'error')
                    return redirect(url_for('.index'))

                # Copy configs
                for config in source_configs:
                    new_config = schema.Config(
                        name=config.name,
                        category=new_category,
                        set_number=new_set_number,
                        value=config.value,
                        param_type=config.param_type,
                        default_value=config.default_value,
                        description=config.description,
                        units=config.units,
                        possible_values=config.possible_values,
                    )
                    self.session.add(new_config)

                self.session.commit()
                flash(f'Configuration set copied to "{new_category}" #{new_set_number}', 'success')
                return redirect(url_for('.edit_set', category=new_category, set_number=new_set_number))

            except Exception as e:
                self.session.rollback()
                flash(f'Error copying configuration: {str(e)}', 'error')

        return self.render('admin/config_set_copy.html',
                           category=category,
                           set_number=set_number)

    @expose('/delete/<category>/<int:set_number>/', methods=['POST'])
    def delete_set(self, category, set_number):
        """Delete a configuration set"""

        try:
            deleted_count = self.session.query(schema.Config).filter(
                schema.Config.category == category,
                schema.Config.set_number == set_number
            ).delete()
            self.session.commit()
            if deleted_count > 0:
                self.session.commit()
                flash(f'Configuration set "{category}" #{set_number} deleted ({deleted_count} parameters)', 'success')
            else:
                flash(f'Configuration set "{category}" #{set_number} not found', 'error')

        except Exception as e:
            self.session.rollback()
            flash(f'Error deleting configuration: {str(e)}', 'error')

        return redirect(url_for('.index'))

    def _create_dynamic_form(self, configs):
        """Create a dynamic WTForms class based on config parameters"""

        class DynamicConfigForm(Form):
            pass

        for config in configs:
            field_name = f"param_{config.ref}"

            # Determine field type based on param_type
            if config.param_type == 'int':
                field = IntegerField(
                    label=config.name,
                    description=config.description,
                    validators=[Optional()]  # Make all fields optional to prevent erasure
                )
            elif config.param_type == 'float':
                field = FloatField(
                    label=config.name,
                    description=config.description,
                    validators=[Optional()]  # Make all fields optional to prevent erasure
                )
            elif config.param_type == 'bool':
                # Check if we should use a dropdown or checkbox
                field = BooleanField(
                    label=config.name,
                    description=config.description
                )
            else:  # string type
                # Check if there are possible values (dropdown)
                if config.possible_values:
                    # Parse possible values and ensure current value is included
                    possible_vals = [v.strip() for v in config.possible_values.split('/') if v.strip()]

                    # Make sure the current value is in the choices (for validation)
                    current_value = config.value.strip() if config.value else ''
                    if current_value and current_value not in possible_vals:
                        possible_vals.append(current_value)

                    # Create choices as (value, label) tuples
                    choices = [(v, v) for v in possible_vals]

                    # Add empty choice for optional fields
                    if not current_value or 'optional' in (config.description or '').lower():
                        choices.insert(0, ('', '-- Select --'))

                    field = SelectField(
                        label=config.name,
                        description=config.description,
                        choices=choices,
                        validators=[Optional()]
                    )
                else:
                    # Long text or regular string
                    if config.description and len(config.description) > 100:
                        field = TextAreaField(
                            label=config.name,
                            description=config.description,
                            validators=[Optional()]
                        )
                    else:
                        field = StringField(
                            label=config.name,
                            description=config.description,
                            validators=[Optional()]
                        )

            # # Add units to description if available
            # if config.units:
            #     current_desc = field.description or ""
            #     field.description = f"{current_desc} (Units: {config.units})".strip()
            #
            # # Add possible values to description if available and not already a select field
            # if config.possible_values and not isinstance(field, SelectField):
            #     current_desc = field.description or ""
            #     field.description = f"{current_desc} Possible values: {config.possible_values}".strip()

            setattr(DynamicConfigForm, field_name, field)

        return DynamicConfigForm
class WorkflowStepView(BaseModelView):
    """Workflow step management"""

    column_list = ['step_id', 'step_name', 'category', 'set_number', 'workflow_id', 'is_active', 'created_at']
    column_searchable_list = ['step_name', 'category', 'workflow_id']
    column_filters = ['category', 'workflow_id', 'is_active']
    column_sortable_list = ['step_id', 'step_name', 'category', 'set_number', 'workflow_id', 'created_at']

    form_columns = ['step_name', 'category', 'set_number', 'workflow_id', 'description', 'is_active']

    column_descriptions = {
        'step_name': 'Unique name for this workflow step',
        'category': 'Configuration category (preprocess, cc, filter, stack, etc.)',
        'set_number': 'Configuration set number',
        'workflow_id': 'Workflow identifier',
        'description': 'Optional description of what this step does',
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'category': lambda view, context, model, name: SafeMarkup(f'<span class="label label-info">{model.category}</span>'),
        'workflow_id': lambda view, context, model, name: SafeMarkup(f'<span class="label label-primary">{model.workflow_id}</span>'),
        'is_active': lambda view, context, model, name: SafeMarkup('✓' if model.is_active else '✗'),
    })

    def __init__(self, session, *args, **kwargs):
        super(WorkflowStepView, self).__init__(schema.WorkflowStep, session, *args, **kwargs)


class WorkflowLinkView(BaseModelView):
    """Workflow link management with smart step selection"""

    column_list = ['link_id', 'from_step', 'to_step', 'link_type', 'is_active', 'created_at']
    column_searchable_list = ['link_type']
    column_filters = ['link_type', 'is_active']
    column_sortable_list = ['link_id', 'from_step_id', 'to_step_id', 'link_type', 'created_at']

    form_columns = ['from_step_id', 'to_step_id', 'link_type', 'is_active']

    column_descriptions = {
        'from_step_id': 'Source workflow step',
        'to_step_id': 'Target workflow step (filtered based on workflow logic)',
        'link_type': 'Type of connection (default, conditional, etc.)',
    }

    # Custom column formatters
    column_formatters = dict(BaseModelView.column_formatters, **{
        'from_step': lambda view, context, model, name: model.from_step.step_name if model.from_step else '',
        'to_step': lambda view, context, model, name: model.to_step.step_name if model.to_step else '',
        'link_type': lambda view, context, model, name: SafeMarkup(f'<span class="label label-default">{model.link_type}</span>'),
        'is_active': lambda view, context, model, name: SafeMarkup('✓' if model.is_active else '✗'),
    })

    def __init__(self, session, *args, **kwargs):
        super(WorkflowLinkView, self).__init__(schema.WorkflowLink, session, *args, **kwargs)

    def create_form(self, obj=None):
        """Override create form to provide smart step selection"""
        form = super().create_form(obj)

        # Get all workflow steps for from_step selection
        from_steps = self.session.query(schema.WorkflowStep).filter(
            schema.WorkflowStep.is_active == True
        ).order_by(schema.WorkflowStep.step_name).all()

        form.from_step_id.choices = [
            (step.step_id, f"{step.step_name} ({step.category}:{step.set_number})")
            for step in from_steps
        ]

        # Initially populate to_step with all steps (will be filtered via JavaScript)
        to_steps = self.session.query(schema.WorkflowStep).filter(
            schema.WorkflowStep.is_active == True
        ).order_by(schema.WorkflowStep.step_name).all()

        form.to_step_id.choices = [
            (step.step_id, f"{step.step_name} ({step.category}:{step.set_number})")
            for step in to_steps
        ]

        return form

    def edit_form(self, obj=None):
        """Override edit form to provide smart step selection"""
        form = super().edit_form(obj)

        # Get all workflow steps
        from_steps = self.session.query(schema.WorkflowStep).filter(
            schema.WorkflowStep.is_active == True
        ).order_by(schema.WorkflowStep.step_name).all()

        form.from_step_id.choices = [
            (step.step_id, f"{step.step_name} ({step.category}:{step.set_number})")
            for step in from_steps
        ]

        to_steps = self.session.query(schema.WorkflowStep).filter(
            schema.WorkflowStep.is_active == True
        ).order_by(schema.WorkflowStep.step_name).all()

        form.to_step_id.choices = [
            (step.step_id, f"{step.step_name} ({step.category}:{step.set_number})")
            for step in to_steps
        ]

        return form

    @expose('/get_valid_targets/<int:from_step_id>')
    def get_valid_targets(self, from_step_id):
        """AJAX endpoint to get valid target steps for a given source step"""
        from_step = self.session.query(schema.WorkflowStep).get(from_step_id)
        if not from_step:
            return jsonify([])

        # Define workflow logic - what steps can follow what
        WORKFLOW_CHAINS = {
            'global': ['preprocess', 'qc'],
            'preprocess': ['cc'],
            'cc': ['filter'],
            'filter': ['stack'],
            'stack': ['mwcs', 'stretching', 'wavelet'],
            'mwcs': ['mwcs_dtt'],
            'stretching': [],  # Terminal step
            'wavelet': ['wavelet_dtt'],
            'wavelet_dtt': [],  # Terminal step
            'mwcs_dtt': [],  # Terminal step
            'qc': []  # Terminal step
        }

        # Get possible next categories for this step's category
        next_categories = WORKFLOW_CHAINS.get(from_step.category, [])

        # Get all steps in those categories
        valid_targets = self.session.query(schema.WorkflowStep).filter(
            schema.WorkflowStep.category.in_(next_categories),
            schema.WorkflowStep.is_active == True,
            schema.WorkflowStep.step_id != from_step_id  # Don't allow self-links
        ).order_by(schema.WorkflowStep.step_name).all()

        return jsonify([
            {
                'id': step.step_id,
                'name': f"{step.step_name} ({step.category}:{step.set_number})"
            }
            for step in valid_targets
        ])

class WorkflowBuilderView(BaseView):
    """
    Visual workflow builder interface
    """

    def __init__(self, session, *args, **kwargs):
        super(WorkflowBuilderView, self).__init__(*args, **kwargs)
        self.session = session

    @expose('/')
    def index(self):
        """Display workflow builder interface"""
        from .api import get_workflow_graph

        # Get available workflows
        workflows = self.session.query(schema.WorkflowStep.workflow_id).distinct().all()
        workflow_ids = [w[0] for w in workflows]

        # Get current workflow (default to 'default')
        current_workflow = request.args.get('workflow', 'default')

        # Get workflow graph data
        graph_data = get_workflow_graph(self.session, current_workflow)

        return self.render('admin/workflow_builder.html',
                           workflows=workflow_ids,
                           current_workflow=current_workflow,
                           graph_data=graph_data)

    @expose('/create_step', methods=['POST'])
    def create_step(self):
        """Create a new workflow step"""
        from .api import create_workflow_step

        step_name = request.form.get('step_name')
        category = request.form.get('category')
        set_number = int(request.form.get('set_number'))
        workflow_id = request.form.get('workflow_id', 'default')
        description = request.form.get('description')

        try:
            step = create_workflow_step(self.session, step_name, category, set_number, workflow_id, description)
            flash(f'Workflow step "{step_name}" created successfully', 'success')
        except Exception as e:
            flash(f'Error creating workflow step: {str(e)}', 'error')

        return redirect(url_for('.index', workflow=workflow_id))

    @expose('/create_link', methods=['POST'])
    def create_link(self):
        """Create a new workflow link"""
        from .api import create_workflow_link

        from_step_id = int(request.form.get('from_step_id'))
        to_step_id = int(request.form.get('to_step_id'))
        link_type = request.form.get('link_type', 'default')

        try:
            link = create_workflow_link(self.session, from_step_id, to_step_id, link_type)
            flash('Workflow link created successfully', 'success')
        except Exception as e:
            flash(f'Error creating workflow link: {str(e)}', 'error')

        return redirect(url_for('.index'))

    @expose('/create_steps_from_configs', methods=['POST'])
    def create_steps_from_configs(self):
        """Create workflow steps automatically from all existing config sets"""
        from .api import create_workflow_steps_from_config_sets

        workflow_id = request.form.get('workflow_id', 'default')

        created_count, existing_count, error_message = create_workflow_steps_from_config_sets(
            self.session, workflow_id
        )

        if error_message:
            flash(f'Error creating workflow steps: {error_message}', 'error')
        else:
            if created_count > 0:
                flash(f'Created {created_count} workflow steps from config sets', 'success')
            if existing_count > 0:
                flash(f'{existing_count} workflow steps already existed', 'info')
            if created_count == 0 and existing_count == 0:
                flash('No config sets found to create workflow steps from', 'warning')

        return redirect(url_for('.index', workflow=workflow_id))

    @expose('/create_links_from_steps', methods=['POST'])
    def create_links_from_steps(self):
        """Create workflow links automatically from existing workflow steps"""
        from .api import create_workflow_links_from_steps

        workflow_id = request.form.get('workflow_id', 'default')

        created_count, existing_count, error_message = create_workflow_links_from_steps(
            self.session, workflow_id
        )

        if error_message:
            flash(f'Error creating workflow links: {error_message}', 'error')
        else:
            if created_count > 0:
                flash(f'Created {created_count} workflow links automatically', 'success')
            if existing_count > 0:
                flash(f'{existing_count} workflow links already existed', 'info')
            if created_count == 0 and existing_count == 0:
                flash('No workflow steps found to create links between', 'warning')

        return redirect(url_for('.index', workflow=workflow_id))

CONFIG_BULK_UPDATE_TEMPLATE = """
{% extends 'admin/base.html' %}

{% block body %}
<div class="container-fluid">
    <h2>Bulk Configuration Update</h2>
    
    <form method="post">
        {% for config in configs %}
        <div class="form-group">
            <label>{{ config.name }} ({{ config.category }})</label>
            <input type="text" name="config_{{ config.ref }}" value="{{ config.value }}" class="form-control">
            <small class="help-block">{{ config.description }}</small>
        </div>
        {% endfor %}
        
        <button type="submit" class="btn btn-primary">Update All</button>
    </form>
</div>
{% endblock %}
"""

MSNOISE_INDEX_TEMPLATE = """
{% extends 'admin/index.html' %}

{% block body %}
<div class="container-fluid">
    <h2>MSNoise Dashboard</h2>
    
    <div class="row">
        <div class="col-md-3">
            <div class="panel panel-primary">
                <div class="panel-heading">Stations</div>
                <div class="panel-body">
                    <h3>{{ stats.get('stations', 0) }}</h3>
                    <p>Active stations</p>
                </div>
            </div>
        </div>
       
        
        <div class="col-md-3">
            <div class="panel panel-success">
                <div class="panel-heading">Workflows</div>
                <div class="panel-body">
                    <h3>{{ stats.get('workflow_steps', 0) }}</h3>
                    <p>Active workflow steps</p>
                </div>
            </div>
        </div>
        
        <div class="col-md-3">
            <div class="panel panel-warning">
                <div class="panel-heading">Jobs</div>
                <div class="panel-body">
                    <h3>{{ stats.get('jobs_todo', 0) }} / {{ stats.get('jobs_total', 0) }}</h3>
                    <p>Todo / Total jobs</p>
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}
"""

import csv
import os
from flask import request, redirect, url_for, flash, render_template_string, render_template
from flask_admin import BaseView, expose
from wtforms import Form, SelectField, StringField, TextAreaField, IntegerField, validators

# ============================================================================
# Main Application Setup
# ============================================================================

def create_admin_app():
    """Create and configure the Flask-Admin application"""

    # Initialize Flask-Admin
    admin = Admin(
        app,
        name='MSNoise Admin',
        template_mode='bootstrap3',
        index_view=MSNoiseAdminIndexView(),
        base_template='admin/base.html'
    )

    # Add model views
    # admin.add_view(FilterView(db_session, name='Filters', category='Configuration'))

    admin.add_view(ConfigSetView(db_session, name='Configuration Sets', endpoint='config_sets', category='Configuration'))
    admin.add_view(StationView(db_session, name='Stations', category='Configuration'))
    admin.add_view(DataAvailabilityView(db_session, name='Data Availability', category='Data'))

    # Add workflow views
    admin.add_view(WorkflowStepView(db_session, name='Workflow Steps', category='Workflow'))
    admin.add_view(WorkflowLinkView(db_session, name='Workflow Links', category='Workflow'))
    admin.add_view(
        WorkflowBuilderView(db_session, name='Workflow Builder', endpoint='workflow_builder', category='Workflow'))

    admin.add_view(JobView(db_session, name='Jobs', category='Processing'))

    admin.add_view(ConfigView(db_session, name='Configuration (Expert)', category='Expert'))
    # Add API view
    # admin.add_view(WorkflowAPIView(name='Workflow API', category='API'))

    # Add custom templates
    app.jinja_env.globals['CONFIG_BULK_UPDATE_TEMPLATE'] = CONFIG_BULK_UPDATE_TEMPLATE
    app.jinja_env.globals['MSNOISE_INDEX_TEMPLATE'] = MSNOISE_INDEX_TEMPLATE


    return admin


# ============================================================================
# Main Entry Point
# ============================================================================

def main(port=5099):
    """Main entry point for the MSNoise admin interface"""

    logger.info("Starting MSNoise Admin Interface")

    # Create admin interface
    admin = create_admin_app()

    # Configure Flask app
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

    # Add template filters
    @app.template_filter('datetime')
    def datetime_filter(value):
        if value is None:
            return ''
        return value.strftime('%Y-%m-%d %H:%M:%S')

    # Error handlers
    @app.errorhandler(404)
    def not_found(error):
        return render_template_string('<h1>Page not found</h1>'), 404

    @app.errorhandler(500)
    def internal_error(error):
        db_session.rollback()
        return render_template_string('<h1>Internal server error</h1>'), 500

    # Cleanup
    @app.teardown_appcontext
    def shutdown_session(exception=None):
        db_session.close()

    # Start the application
    try:
        app.run(host='0.0.0.0', port=port, debug=True)
    except KeyboardInterrupt:
        logger.info("Shutting down MSNoise Admin Interface")
    finally:
        db_session.close()


if __name__ == '__main__':
    main()
